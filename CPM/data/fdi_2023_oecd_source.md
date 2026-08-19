# fdi_2023_oecd.csv の出典

## 書誌情報

OECD. *Foreign direct investment positions by partner country*.
OECD Data Explorer / SDMX. https://sdmx.oecd.org/ （2026年8月19日取得）

## 取得条件

| 項目 | 値 |
|---|---|
| データフロー | `OECD.DAF.INV,DSD_FDI@DF_FDI_POS_CTRY,1.0` |
| 年 | 2023年末 |
| MEASURE | `LE_FA_F`（FDI positions - total） |
| MEASURE_PRINCIPLE | **`DO`（Directional principle: outward）** |
| ACCOUNTING_ENTRY | `NET_FDI` |
| TYPE_ENTITY | `ALL`（全事業体） |
| LEVEL_COUNTERPART | `IMC`（即時相手国） |
| SECTOR / ACTIVITY | `S1` / `_T` |
| 単位 | 米ドル（`UNIT_MULT` を乗じて実額に換算済み） |
| 対象国 | USA, JPN, KOR, DEU, FRA, NLD, BEL, LUX |

生レスポンスの要約は `fdi_2023_raw.json`（取得時刻を含む）。

## 行と列の向き（確定済み）

**行が投資国、列が投資先である。** すなわち `T[i][j]` は国 i から国 j への
対外直接投資ポジションであり、各セルは投資国自身の報告値である。

この向きは OECD のコードリストと実値の両方で確認した。

- コードリスト：`DO` = *Directional principle: outward*、`DI` = *Directional principle: inward*
- 実値による照合（2023年、十億ドル）

  | 系列 | REF_AREA | COUNTERPART | 値 | 公表値との整合 |
  |---|---|---|---|---|
  | `DO` | JPN | USA | 705.7 | 日本の対米投資（約700）と一致 |
  | `DO` | USA | JPN | 59.9 | 米国の対日投資 |
  | `DI` | USA | JPN | 653.5 | 上の JPN→USA を米国側が報告したもの |

当初 `DI`（内向き）で取得していたため行と列が逆になっていたが、`DO` に修正した。
これにより、貿易データを用いた別稿（行＝輸出国、送り手側の報告で統一）と
同じ方針、すなわち「送り手側が報告した値で全セルを揃える」が保たれる。

## 年の選択

`DO` 系列は2020年から2023年まで、8か国56方向すべてが欠損なく、
かつすべて正の値である。最新である2023年を採用した。
2020〜2022年も同条件を満たすため、頑健性の確認に使える。

なお `DI`（内向き）系列は2023年の即時相手国ベースがほぼ未公表で、
双方向がそろうペアが存在しない。`DO` を用いる理由の一つでもある。

## 検算

`fdi_2023_oecd.csv` の値は、本プロジェクトの初期資料に含まれていた
FDI 行列（USA 行：63, 36, 193, 101, 980, 67, 532）とほぼ一致する
（本データ：59.9, 34.4, 201.3, 101.4, 995.4, 64.4, 522.4、単位は十億ドル）。
