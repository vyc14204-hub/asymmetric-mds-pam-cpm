# ============================================================================
#  PAM（極座標表現）と CPM（循環成分）による非対称MDS
#
#  上から下へ順に読める形で書いてある。関数にまとめず、似た処理は書き下してある。
#  冗長だが、どの行が何をしているかを追えることを優先している。
#
#  前提：R Canvas のシート fdi_2023_oecd
#        （OECD 対外直接投資ポジション 2023年末、行=投資国、列=投資先、単位=米ドル）
#  必要：install.packages(c("ggplot2", "ggrepel", "smacof"))
# ============================================================================

library(ggplot2)
library(ggrepel)
library(smacof)


# ============================================================================
#  1. データを行列にする
# ============================================================================

fdi <- fdi_2023_oecd
rownames(fdi) <- fdi$country
T <- data.matrix(fdi[, -1])     # 1列目（国名）を除く

diag(T) <- NA                   # 自国に対する値は無い

nm <- rownames(T)
n  <- nrow(T)                   # 8

print(T / 1e9)                  # 10億ドル単位で表示


# ============================================================================
#  2. 非類似度行列を2通り作る
# ============================================================================
# PAM と CPM で必要な変換が違うので、別々に作る。
#
#   PAM 用：d = 100 / T（T は10億ドル単位）
#       偏角 arctan(a/s) は s の符号に依存するので、非類似度が正である必要がある。
#       対数変換だと s が負になり偏角の意味が壊れるため、逆数変換を使う。
#
#   CPM 用：d = -log(T)
#       歪対称成分 A は非類似度への加法定数に不変なので、対数変換でよい。
#       循環成分 C は A の線形変換なので同じ不変性をもつ。

Tb <- T / 1e9                   # 10億ドル単位

D_rec <- 100 / Tb               # PAM 用
diag(D_rec) <- 0

D_log <- -log(T)                # CPM 用
diag(D_log) <- 0


# ============================================================================
#  3. PAM：対称成分と歪対称成分から偏角を作る
# ============================================================================
# 対称成分 s と歪対称成分 a を実部・虚部とみなした複素数 z = s + i a の偏角。
#     θ_ij = arctan(a_ij / s_ij)
# atan2 を使うと符号と象限を正しく扱える。

S_rec <- (D_rec + t(D_rec)) / 2
A_rec <- (D_rec - t(D_rec)) / 2

theta <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      theta[i, j] <- atan2(A_rec[i, j], S_rec[i, j])
    }
  }
}

# 歪対称になっているか確認（θ_ij = -θ_ji のはず）
cat("θ は歪対称か :", all(abs(theta + t(theta)) < 1e-12), "\n")

# 度に直して大きい順に見る
theta_deg <- theta * 180 / pi
cat("\n偏角が大きいペア（度）\n")
for (i in 1:(n-1)) for (j in (i+1):n) {
  if (abs(theta_deg[i, j]) > 35) {
    cat(sprintf("   %s-%s  %.2f\n", nm[i], nm[j], theta_deg[i, j]))
  }
}

# θ^2 を二乗非類似度とみなすので、非類似度そのものは |θ|
PAM <- abs(theta)


# ============================================================================
#  4. PAM の布置
# ============================================================================
# まず古典的MDSを試して、負の固有値がどれくらい出るかを見る。

fit_pam_c <- cmdscale(as.dist(PAM), k = 2, eig = TRUE)
ev <- fit_pam_c$eig
cat("\nPAM 古典的MDS：2次元説明率",
    round(sum(ev[1:2]) / sum(abs(ev)) * 100, 1), "%\n")
cat("PAM 最小固有値 / 最大固有値 :",
    round(-min(ev) / max(ev) * 100, 1), "%\n")

# 負の固有値が大きいので、SMACOF を本採用する
set.seed(123)
fit_pam <- mds(as.dist(PAM), ndim = 2, type = "ratio")
X_pam <- fit_pam$conf

cat("PAM SMACOF stress-1 :", round(fit_pam$stress, 4), "\n")


# ============================================================================
#  5. CPM：歪対称成分を勾配と循環に分ける
# ============================================================================
# 対数変換の側で歪対称成分を作る。

A <- (D_log - t(D_log)) / 2

# 勾配 φ_i：その対象が他の対象に対してもつ歪対称成分の平均
phi <- rowSums(A) / n

cat("\n勾配 φ\n")
print(round(sort(phi), 3))

# 勾配だけで説明される非対称量 G_ij = φ_i - φ_j
G <- outer(phi, phi, "-")

# 残りが循環成分
C <- A - G
diag(C) <- 0

# C の性質を確認する
cat("\nC は歪対称か   :", all(abs(C + t(C)) < 1e-10), "\n")
cat("C の行和はゼロか :", all(abs(rowSums(C)) < 1e-10), "\n")

# 平方和の分け前（A のうちどれだけが勾配で、どれだけが循環か）
ss_A <- sum(A^2); ss_G <- sum(G^2); ss_C <- sum(C^2)
cat(sprintf("平方和  A = %.3f = G %.3f (%.1f%%) + C %.3f (%.1f%%)\n",
            ss_A, ss_G, ss_G/ss_A*100, ss_C, ss_C/ss_A*100))


# ============================================================================
#  6. CPM を計算する
# ============================================================================
#     d_CPM^2(i,j) = sum_{k=1}^{n} (c_ki - c_kj)^2
#
# C の第 i 列と第 j 列のユークリッド距離。APM を A ではなく C に対して作ったもの。

CPM2 <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      CPM2[i, j] <- sum((C[, i] - C[, j])^2)
    }
  }
}

CPM <- sqrt(CPM2)


# ============================================================================
#  7. 分解の検算
# ============================================================================
# C も歪対称かつ対角ゼロなので、A のときと同じ分解が成り立つはず。
#     d_CPM^2 = 2 c_ij^2 + sum_{k != i,j} (c_ki - c_kj)^2

CDD2 <- C^2                                    # 直接項
CTP2 <- matrix(0, n, n, dimnames = list(nm, nm))   # 第三者項

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      k <- setdiff(1:n, c(i, j))
      CTP2[i, j] <- sum((C[k, i] - C[k, j])^2)
    }
  }
}

cat("C についても分解が成立するか :",
    all(abs(CPM2 - (2 * CDD2 + CTP2)) < 1e-10), "\n")

share_c <- 2 * CDD2 / CPM2
diag(share_c) <- NA
cat(sprintf("直接項の寄与 : 最小 %.1f%%  最大 %.1f%%  平均 %.1f%%\n",
            min(share_c, na.rm = TRUE), max(share_c, na.rm = TRUE),
            mean(share_c, na.rm = TRUE)) )


# ============================================================================
#  8. CPM の布置
# ============================================================================
# CPM は C の列ベクトル間のユークリッド距離なので、必ずユークリッド距離行列。
# したがって古典的MDSで問題ない。

fit_cpm <- cmdscale(as.dist(CPM), k = 2, eig = TRUE)
X_cpm <- fit_cpm$points

ev2 <- fit_cpm$eig
cat("CPM 2次元説明率 :", round(sum(ev2[1:2]) / sum(abs(ev2)) * 100, 1), "%\n")
cat("CPM 最小固有値  :", round(min(ev2), 6), " （0以上ならユークリッド）\n")


# ============================================================================
#  9. 作図
# ============================================================================
# 2つとも同じ作りなので2回書く。coord_equal() は MDS 図では必須。

# --- 図：PAM-MDS ---
df_pam <- data.frame(country = nm, Dim1 = X_pam[, 1], Dim2 = X_pam[, 2])

p_pam <- ggplot(df_pam, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "PAM-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p_pam)

# --- 図：CPM-MDS ---
df_cpm <- data.frame(country = nm, Dim1 = X_cpm[, 1], Dim2 = X_cpm[, 2])

p_cpm <- ggplot(df_cpm, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "CPM-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p_cpm)


# ============================================================================
#  10. 表
# ============================================================================
# 28ペア（8国から2つ選ぶ組み合わせ）を1行ずつ並べる。

# --- PAM の表：対称成分、歪対称成分、偏角 ---
tbl_pam <- data.frame()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    tbl_pam <- rbind(tbl_pam, data.frame(
      pair      = paste(nm[i], nm[j], sep = "-"),
      s         = round(S_rec[i, j], 4),
      a         = round(A_rec[i, j], 4),
      theta_deg = round(theta_deg[i, j], 2)
    ))
  }
}

print(tbl_pam[order(-abs(tbl_pam$theta_deg)), ], row.names = FALSE)

# --- CPM の表：循環成分、第三者項、CPM、直接項の寄与 ---
tbl_cpm <- data.frame()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    tbl_cpm <- rbind(tbl_cpm, data.frame(
      pair   = paste(nm[i], nm[j], sep = "-"),
      c_abs  = round(abs(C[i, j]), 3),
      third  = round(sqrt(CTP2[i, j]), 3),
      CPM    = round(CPM[i, j], 3),
      share  = round(2 * CDD2[i, j] / CPM2[i, j] * 100, 1)
    ))
  }
}

print(tbl_cpm[order(-tbl_cpm$share), ], row.names = FALSE)
