# 内容
本コードは、銀河のイメージ画像をシミュレーションし、楕円率を測るコードです。

宇宙に広がるダークマターや銀河による重力レンズ効果で遠方銀河の形が変形する現象(弱い重力レンズ効果)は、ダークマターを含む宇宙の全物質分布を調べることができる貴重な観測手法です。
重力レンズ効果の測定には非常に多くの銀河の形を精確に測る必要があります。
しかし、銀河の種類によって輝度プロファイルは大きく異なるうえ、地上から観測した銀河の像は大気ゆらぎや大気分散による像の広がり(PSF)やノイズの影響を受けており、
精密かつ高速に銀河の形を測るのは容易ではありません。

本コードはSpergel 2010の方法を使い銀河の形(楕円率)を測定します。
この方法では、銀河の輝度分布を第3種ベッセル関数でフィットし、PSFの影響を補正することで銀河の形(楕円率)を測定するコードです。
第3種ベッセル関数は、楕円銀河(ドゥ・ボークルール則)や渦巻銀河のディスク(指数則)の輝度プロファイルをよくフィットできるうえ、
フーリエ成分の解析的な表式があるため高速な処理ができます。

# Refererences

Analytical Galaxy Profiles for Photometric and Lensing Analysis  
D. N. Spergel  
The Astrophysical Journal Supplement, Volume 191, Issue 1, pp. 58-65 (2010)
