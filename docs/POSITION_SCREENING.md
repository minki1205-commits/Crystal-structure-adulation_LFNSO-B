# XRD 그림 / 추가 TIFF를 이용한 위치 가설 비교

## 입력과 방법

- 사용자가 공유한 A/B XRD 그림을 눈으로 판독했다. B의 주요 d-spacing은 약 4.89, 4.50, 4.03, 3.62, 3.08, 2.57 Å로 기록했다. 보이는 A 분리 피크는 약 4.30, 3.85, 3.31, 2.84, 2.56 Å이다. 겹친 피크를 A의 독립 관측으로 추가하지 않았다.
- 수동 판독 범위는 ±0.02 Å로 설정했다. 통계적 표준오차가 아니고 기기 오차도 포함하지 않는다. Inline 그림은 실행 환경에 raster 파일로 제공되지 않아 자동 색상 digitization을 했다고 주장하지 않는다. XRD raw profile/파장은 아직 없다. d-spacing 비교에는 파장이 필요하지 않지만 강도 피팅·2θ 기반 분석·정밀 격자 refinement는 하지 않았다.
- 그림의 q≈2.2, 2.48 nm⁻¹는 우선 g=1/d라는 절대 반사 위치로 취급했다. 다른 hkl의 |g| 차이를 modulation vector로 간주하지 않았다.
- TIFF 두 개는 각각 696² / 683²의 8-bit palette display export이다. 10 nm⁻¹ scale bar 길이 126/148 pixel로 0.079365/0.067568 nm⁻¹/pixel을 얻었다. 양 끝 합계 2 pixel의 판독 범위는 약 1.59%/1.35% scale 범위이며 실제 calibration uncertainty를 대신하지 않는다. Scale bar·문자 영역은 ROI에서 제외했다. 표시용 밝기는 raw 전자 counts가 아니다.

## 1. 고정 시작 격자의 XRD 가설 비교

격자 a=5.1706, b=8.9382, c=5.1635 Å, β=109.49°를 고정했다. B의 세 피크를 각각 parent (110), (-111), (021)의 +δb* satellite로 대응시킨 조건부 비교다. 이 대응은 검증할 가설이며 독립적으로 확정된 indexing이 아니다.

| B 그림 d / Å | δ=1/4 예측 / Å | 차이 / Å |
|---|---|---|
| 4.03 | 4.02758 | -0.00242 |
| 3.62 | 3.63389 | +0.01389 |
| 3.08 | 3.07769 | -0.00231 |

- 위치 RMS: δ=1/4 0.00825 Å, δ=1/3 0.07060 Å, δ=1/6 0.07546 Å.
- 자유 δ의 최소값은 약 0.25266이나, 그림과 starting metric 수준에 비해 이 소수 자릿수는 정밀도를 의미하지 않는다.
- 이 비교는 원자 구조의 순위가 아니다. Occupancy/stacking/displacement의 강도와 소멸 조건은 계산하지 않았다. 12b supercell의 3차 harmonic 등 동일 q를 만드는 다른 구조도 배제되지 않는다.

## 2. 격자 변화와 q 사이의 모호성

시작 격자의 모든 길이를 각각 0.99, 1.00, 1.01배로 바꾸는 진단에서 최적 δ는 약 0.211, 0.253, 0.293으로 바뀐다. ±1%는 측정 오차로 가정한 값이 아니라 민감도 확인용 범위다. 따라서 그림만으로 정확히 4b라고 확정할 수 없다.

B의 4.89/4.50 Å 피크가 각각 (001)/(020)이라는 조건에서 monoclinic unique-b metric은

`1/d(021)^2 = 1/d(001)^2 + 1/d(020)^2`

를 만족해야 한다. 예상 d(021)는 3.311 Å (수동 판독범위를 전파하면 약 3.297–3.325 Å)로, B의 3.08 Å 피크와 다르다. 따라서 이 세 피크의 동일-hkl 대응을 모두 유지한 채 단순 metric 변화만으로 설명하는 해석은 맞지 않는다. 다른 hkl, superstructure, 상/도메인 중첩은 별도 검토해야 한다.

## 3. DM3 회절의 경험적 q 추정

DP01_101.dm3의 강한 두 반사에서 경험적 reciprocal basis u,v를 만들고 C-centering 형태의 정수 parity로 22개 parent-like spot을 정합했다. 아직 결정학적으로 확정된 hkl basis는 아니다.

- |u|≈3.14397, |v|≈1.07059 nm⁻¹, 두 축 각도≈89.921°.
- 선택 조건을 만족하는 satellite-like spot 22개의 중앙 잔차 δ_v≈0.28823; robust scatter≈0.00653. 반대 spot 및 공통 basis 때문에 이들을 독립 표본으로 세거나 scatter를 최종 오차로 해석하지 않는다.
- transverse residual δ_u≈-0.00731, 영상 좌표 q_yx≈(0.04768,-0.30570) nm⁻¹, |q|≈0.3094 nm⁻¹.
- parent |k|<4 그룹 중앙값≈0.286, 4≤|k|<8 그룹≈0.291. 위치/왜곡/선택에 따른 영향을 추가 확인해야 한다.
- 이는 v≈b*라는 조건하에서 1/4과 1/3 사이에 있다. 2/7, 3/10 같은 rational approximant나 비정합 변조를 곧바로 결정하지 않는다. 원본에서 검증된 3D q가 아니다.

## 4. XRD와 ED를 함께 설명할 수 있는가?

경험적 ED의 δ≈0.28823을 고정하고, XRD 세 피크에 대해 전체 길이 scale 하나만 피팅하면 scale≈1.008416 (+0.842%)에서 예측 d≈4.0220, 3.6354, 3.0722 Å, RMS≈0.0110 Å이다. 피팅에 사용하지 않은 (001)/(020)의 예상 d≈4.9086/4.5067 Å도 그림 판독범위 부근이다. 그러므로 XRD와 ED의 상대 q가 반드시 모순이라고 할 수 없다.

그러나 이 조건의 metric b*≈1.1095 nm⁻¹와 ED에서 얻은 |v|≈1.0706 nm⁻¹에는 약 3.6% scale 차이가 남는다. 이를 metadata calibration에 몰래 보정하지 않았다. Calibration, 시료 내 metric 차이, indexing/setting 문제를 분리해야 한다. 이 피팅은 완전한 joint refinement가 아니다.

파일명 [101]의 축 표기도 확인할 필요가 있다. 측정 |u|/|v|≈2.937은 현재 unique-b convention에서 [101]의 예상 g(10-1)/b*≈2.118보다 [10-1]의 g(101)/b*≈2.997에 가깝다. 이는 v=b*, u가 해당 기본 벡터라는 조건부 점검이다. 단위격자 축 부호나 setting 차이일 수 있으므로 파일명/zone을 임의로 바꾸지 않았다.

## 5. TIFF 결과와 남은 범위

- [001] TIFF: 139개 자동 피크 후보, 반대편 반사쌍 49개로 중심 추정.
- [110] TIFF: 160개 자동 피크 후보(설정한 상한에 도달), 반대편 반사쌍 21개로 중심 추정. 상한·threshold 의존성 때문에 전체 반사 개수가 아니다.
- 미리보기에서 좌표·center를 확인했고 scale bar를 spot으로 세지 않았다. 축소/밝기 포화로 작은 satellite의 분리 및 intensity 정량은 제한된다.
- TIFF 위치 자료는 기존 이상값을 가진 DM4의 raw intensity를 대체하지 않는다. 두 자료의 provenance를 구분했다.
- LFNSO 화학적 parent CIF/기존 모델 코드가 아직 없어 구조 생성·강도 비교·STEM 단계로 진입하지 않았다. XRD/ED gate는 not evaluated이다.

## 재현과 출력

저장소 루트에서:

```bash
MPLCONFIGDIR=/workspace/lfnso-cache/matplotlib PYTHONPATH=src /workspace/lfnso-venv/bin/python -m lfnso.position_screening configs/position_screening.json
```

이미 cycle 디렉터리가 존재하면 output에 새 이름을 지정한다. 원본/기존 결과를 덮어쓰지 않는다. 현재 결과는 `03_structure_model_search/position_screening_02/`이다. 설정에 수동 XRD 판독값·선택한 reciprocal-basis seed·TIFF ROI와 scale bar 좌표를 모두 저장했다.

다음에 받으면 가장 도움이 되는 것은 기존 LFNSO CIF/생성 코드와 hkl 또는 축 방향이 표시된 회절 indexing 그림이다. XRD raw 자료가 도착하면 수동 판독값을 교체하고 격자–q 모호성을 다시 검증한다.
