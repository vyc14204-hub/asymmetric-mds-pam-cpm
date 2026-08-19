# fdi_2021_oecd.csv の出典

## 書誌情報

OECD. *Foreign direct investment positions by partner country*.
OECD Data Explorer / SDMX. https://sdmx.oecd.org/ （2026年8月19日取得）

## 取得条件

| 項目 | 値 |
|---|---|
| エンドポイント | `https://sdmx.oecd.org/public/rest/data/` |
| データフロー | `OECD.DAF.INV,DSD_FDI@DF_FDI_POS_CTRY,1.0` |
| 年 | 2021年末 |
| MEASURE | `LE_FA_F`（ポジション） |
| MEASURE_PRINCIPLE | `DI`（方向性原則） |
| ACCOUNTING_ENTRY | `NET_FDI` |
| TYPE_ENTITY | `ALL`（全事業体） |
| LEVEL_COUNTERPART | `IMC`（即時相手国） |
| SECTOR / ACTIVITY | `S1` / `_T` |
| 単位 | 米ドル（UNIT_MULT=6 を乗じて実額に換算済み） |
| 対象国 | USA, JPN, KOR, DEU, FRA, NLD, BEL, LUX |

生レスポンスの要約は `fdi_2021_raw.json`（取得時刻を含む）。

## 年を2021年にした理由

- **2023年**：公表が最終投資国（ULT）ベースのみで、即時相手国ベースがほぼ空。
  8か国について双方向がそろうペアが **0組** であり、非対称行列が構成できない。
- **2022年**：28ペアすべて双方向がそろうが、JPN→BEL が **−0.03（10億ドル）** と負値。
  逆投資により FDI ポジションは負になりうるが、`d = 100/T` や `-log T` が定義できない。
- **2021年**：28ペアすべて双方向がそろい、対角以外すべて正。**これを採用した。**

2020年も同条件を満たすので、頑健性の確認に使える。

## 注意：行と列の向き

`REF_AREA` を行、`COUNTERPART_AREA` を列としているが、方向性原則・`NET_FDI` の
組み合わせにおいて「行の国が列の国へ投資した額」と読めるかは未確認である。
たとえば USA 行の NLD 列（705.5、10億ドル）や JPN 列（653.5）は、
米国から見た対外投資というより、両国から米国への投資規模に近い水準である。
**本文で投資の向きに言及する前に、OECD の系列定義を確認すること。**
数学的な処理（歪対称成分の抽出、勾配・循環分解、偏角）は向きの解釈に依存しないため、
分析結果そのものは影響を受けない。
