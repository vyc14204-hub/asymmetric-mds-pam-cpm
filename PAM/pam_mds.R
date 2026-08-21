# ============================================================================
#  PAM（極座標表現）による非対称MDS
#
#  上から下へ順に読める形で書いてある。関数にまとめていない。
#  冗長だが、どの行が何をしているかを追えることを優先している。
#
#  前提：data/fdi_2023_oecd.csv（このリポジトリに同梱）
#        （OECD 対外直接投資ポジション 2023年末、行=投資国、列=投資先、単位=10億米ドル）
#  必要：install.packages(c("ggplot2", "ggrepel", "smacof"))
# ============================================================================

library(ggplot2)
library(ggrepel)
library(smacof)


# ============================================================================
#  1. データを行列にする
# ============================================================================

fdi <- read.csv("data/fdi_2023_oecd.csv", fileEncoding = "UTF-8-BOM")
rownames(fdi) <- fdi$country
T <- data.matrix(fdi[, -1])     # 1列目（国名）を除く

diag(T) <- NA                   # 自国に対する値は無い

nm <- rownames(T)
n  <- nrow(T)                   # 8

print(T)                        # 単位は10億米ドル


# ============================================================================
#  2. 投資額を非類似度に変換する
# ============================================================================
#     δ_ij = 100 / T_ij   （T_ij は10億米ドル単位）
#
# 偏角 arctan(a/s) は s の符号に依存するので、非類似度が正である必要がある。
# 対数変換だと s が負になり偏角の意味が壊れるため、逆数変換を使う。

Delta <- 100 / T
diag(Delta) <- 0


# ============================================================================
#  3. 対称成分と歪対称成分に分ける
# ============================================================================

S <- (Delta + t(Delta)) / 2
A <- (Delta - t(Delta)) / 2


# ============================================================================
#  4. 偏角を計算する
# ============================================================================
# 対称成分 s を実部、歪対称成分 a を虚部とみなした複素数 z = s + i a の偏角。
#     θ_ij = arctan(a_ij / s_ij)
# atan2 を使うと符号と象限を正しく扱える。

theta <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      theta[i, j] <- atan2(A[i, j], S[i, j])
    }
  }
}

# 歪対称になっているか確認（θ_ij = -θ_ji のはず）
cat("θ は歪対称か :", all(abs(theta + t(theta)) < 1e-12), "\n")

theta_deg <- theta * 180 / pi

cat("\n偏角が大きいペア（度）\n")
for (i in 1:(n-1)) for (j in (i+1):n) {
  if (abs(theta_deg[i, j]) > 30) {
    cat(sprintf("   %s-%s  %.2f\n", nm[i], nm[j], theta_deg[i, j]))
  }
}

# 偏角の絶対値 |θ| を非類似度とする
PAM <- abs(theta)


# ============================================================================
#  5. 布置を求める
# ============================================================================
# まず古典的MDSを試して、負の固有値がどれくらい出るかを見る。

fit_c <- cmdscale(as.dist(PAM), k = 2, eig = TRUE)
ev <- fit_c$eig

cat("\n古典的MDS 2次元説明率 :", round(sum(ev[1:2]) / sum(abs(ev)) * 100, 1), "%\n")
cat("最小固有値 / 最大固有値 :", round(-min(ev) / max(ev) * 100, 1), "%\n")

# 負の固有値が大きいので SMACOF を本採用する
set.seed(123)
fit <- mds(as.dist(PAM), ndim = 2, type = "ratio")
X <- fit$conf

cat("SMACOF stress-1 :", round(fit$stress, 4), "\n")


# ============================================================================
#  6. 作図
# ============================================================================
# coord_equal() は MDS 図では必須。縦横の縮尺が違うと距離が正しく見えない。

df <- data.frame(country = nm, Dim1 = X[, 1], Dim2 = X[, 2])

p <- ggplot(df, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "PAM-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p)


# ============================================================================
#  7. 表3のもとになる一覧
# ============================================================================
# 28ペア（8国から2つ選ぶ組み合わせ）を1行ずつ並べる。

tbl <- data.frame()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    tbl <- rbind(tbl, data.frame(
      pair      = paste(nm[i], nm[j], sep = "-"),
      s         = round(S[i, j], 4),
      a         = round(A[i, j], 4),
      theta_deg = round(theta_deg[i, j], 2)
    ))
  }
}

print(tbl[order(-abs(tbl$theta_deg)), ], row.names = FALSE)


# ============================================================================
#  8. 布置の解釈に使う量
# ============================================================================
# 各対象が関与する偏角の総和。第1次元が何に対応しているかを見るため。

total <- colSums(abs(theta_deg))

cat("\n国別 |θ| の合計（度）\n")
print(round(sort(total, decreasing = TRUE), 1))

cat("第1次元との相関 :", round(cor(total, X[, 1]), 3), "\n")


# ============================================================================
#  8. 検証：偏角の上限は 45 度である（90 度ではない）
# ============================================================================
# 非類似度が非負なら |a_ij| = |δ_ij - δ_ji| / 2 <= (δ_ij + δ_ji) / 2 = s_ij
# が恒等的に成り立つ。したがって |θ| は 45 度を超えない。

off <- !diag(n)
cat("\n--- 偏角の値域 ---\n")
cat("  |a| <= s がすべて成り立つか :", all(abs(A[off]) <= S[off] + 1e-12), "\n")
cat("  観測された最大 |θ|          :", round(max(abs(theta_deg)), 3), "度\n")
cat("  理論上の上限                : 45 度\n")


# ============================================================================
#  9. 検証：この変換のもとで PAM は往復比のみの関数である
# ============================================================================
#     δ_ij = c / T_ij を代入すると a/s = (T_ji - T_ij) / (T_ji + T_ij) となり、
#     θ_ij = arctan(T_ji / T_ij) - 45度
# が厳密に成り立つ。定数 c は約分されて消え、投資規模は偏角に影響しない。

ident <- matrix(0, n, n, dimnames = list(nm, nm))
for (i in 1:n) for (j in 1:n) if (i != j)
  ident[i, j] <- atan(T[j, i] / T[i, j]) * 180 / pi - 45

cat("\n--- 恒等式 θ_ij = arctan(T_ji / T_ij) - 45度 ---\n")
cat("  最大絶対誤差 :", format(max(abs(theta_deg[off] - ident[off])), digits = 3), "\n")


# ============================================================================
#  10. 検証：|a| を使う方法とは重視するペアが異なる
# ============================================================================

u   <- upper.tri(theta)
cmp <- data.frame(pair       = outer(nm, nm, paste, sep = "-")[u],
                  theta_deg  = round(abs(theta_deg[u]), 2),
                  rank_theta = rank(-abs(theta_deg[u])),
                  abs_a      = round(abs(A[u]), 3),
                  rank_a     = rank(-abs(A[u])))

cat("\n--- |θ| 上位5ペア ---\n")
print(head(cmp[order(cmp$rank_theta), ], 5), row.names = FALSE)
cat("\n--- |a| 上位5ペア ---\n")
print(head(cmp[order(cmp$rank_a), ], 5), row.names = FALSE)
cat("\n順位相関（Spearman）:", round(cor(cmp$rank_theta, cmp$rank_a), 3), "\n")
