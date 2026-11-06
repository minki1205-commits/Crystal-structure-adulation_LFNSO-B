# 실제 ED–후보 위치 매칭 및 공동 검토

2026-11-06 업데이트. 기존 경험적 q/기하 반사표 단계에서 실제 검출 좌표와 후보 계산 패턴의 정합으로 확장했다. 최신 occupational_search_03의 **후보 12개 + ordered control 2개**, 영상 네 장, zone 대안까지 **70개 조합**을 비교했다. 이전 탐색 종점 8개는 이 비교에 다시 넣지 않았다.

## 현재 process

1. 입력 정리 및 추출: 접근 가능한 XRD/ED/STEM 자료를 보존하고 peak 후보, metadata, hash를 기록했다. 일부 기기 조건과 이전 Windows 작업본은 여전히 없다.
2. 가설 생성: 앞선 32개 모델 및 조성 보존 자리 교환 탐색 종점 20개를 생성했다. 구조/산화수의 확정이나 에너지 relaxation은 아니다.
3. XRD screening: 실험 프로파일과 비교했고 11b 후보 일부가 지정 superstructure 구간에서 우세했다. 불필요 반사와 강도 불일치로 gate 통과는 성립하지 않았다.
4. **ED 위치 정합: 이번에 구현·실행했다.** 회전/배율/중심 정합, 후보별 일대일 매칭, 위치 허용오차·반사 threshold 민감도, 회전 대조군, 추가 방향별 왜곡 진단을 포함한다.
5. XRD+ED 공동 구조 선택: 미완료. 점수 정의와 불확실성이 달라 임의 가중합으로 최종 1위를 정하지 않았다.
6. STEM multislice: 미실행. 단순 원자 투영이나 Z-weight proxy를 STEM 시뮬레이션으로 부르지 않는다.

## 입력 및 좌표 처리

- `DP01_101.dm3`, `dp3_110.dm3`: 원본 검출 좌표 및 metadata scale 0.016421644 nm⁻¹/pixel.
- `[001] OneView` 및 `[110] OneView` TIFF: 표시용 export의 검출 좌표, scale bar 환산값 10/126 및 10/148 nm⁻¹/pixel. 정량 electron counts로 취급하지 않는다.
- 기존 detection ROI와 0.6–12 nm⁻¹ 반경의 교집합을 사용했다. Beam 중심 제외와 2 pixel 이내 중복 centroid 제거 후 각각 61/87/118/120개 좌표가 남았다. 모든 점이 확정된 Bragg reflection이라는 뜻은 아니다.
- 원본의 peak extraction을 바꾼 것은 아니다. 각 match에 원래 peak index를 보존했다. 저장 FFT는 원본 HAADF와 중복이므로 독립 ED 영상으로 추가하지 않았다. 이상값이 있는 DM4는 사용하지 않았다.

## 계산과 정합

각 후보의 정확한 parent→supercell 변환을 사용했다. Zone은 parent direct basis로 주어지고 `h*u + (k_super/n)*v + l*w = 0`인 strict zero-order Laue-zone(ZOLZ)만 생성한다. Reciprocal convention은 2π를 제외한 nm⁻¹이다.

전자 산란인자는 neutral-atom Mott–Bethe/Born 근사이며 X-ray 산란 계수로부터 계산했다. Fe/Ni의 .67/.33 점유율을 구조 인자 합에 명시적으로 곱했다. 이 점은 라이브러리 `TEMCalculator.cell_scattering_factors`가 site occupancy를 사용하지 않는 것과 구분한다. 원소별 계수는 라이브러리 계산과 교차 검증했다. 두께, tilt, Ewald excitation, 다중 산란, 열진동, 이온 산란인자를 정련하지 않았다.

기본 정합은 영상·metric·zone마다 한 번 수행한다. 처음 32개 검출 좌표의 65% trimmed positional loss를 parent C-centered reciprocal grid에 맞추고, 회전·mirror·등방 배율·중심 이동을 탐색한다. 배율 범위 0.9–1.1, 중심 ±3 pixel은 **진단 범위이지 실측 오차 범위가 아니다**. 같은 metric/zone의 후보는 같은 pose를 사용하므로 후보마다 자유롭게 영상을 왜곡하지 않는다. [101] 파일은 [101]/[10-1] 가설을 모두 남겼다.

후보별 계산 반사는 비교 ROI 내 최고 kinematic intensity의 0.1/1/5% threshold로 선택했다. 이는 experimental detection limit이 아니다. 위치 허용오차는 DM3 0.05/0.08/0.12 nm⁻¹, TIFF 0.16/0.24/0.32 nm⁻¹를 시험했다. 허용오차 이내 일대일 대응 수를 우선 최대화한 뒤 거리를 최소화한다. 하나의 반사를 여러 실험 점에 중복 대응시키지 않는다.

Recall=matched/observed, precision=matched/predicted, F1=2·matched/(observed+predicted)를 기록했다. Dense supercell이 우연히 많은 점을 설명하는 효과를 줄이기 위해 precision과 동일 후보를 13/29/47/71° 추가 회전한 대조 결과도 기록했다. 이는 통계적 p-value 검정이 아니다. 미검출 spot은 threshold/포화/beam stop/동역학 영향이 있으므로 단순히 실제 소멸로 판정하지 않는다.

## 결과 및 제한

조건부 B metric의 trimmed registration RMS는 [101] 가설 0.419 nm⁻¹, [10-1] 가설 0.126 nm⁻¹이다. 현재 setting에서는 후자가 더 잘 맞지만 파일명이나 사용자 zone을 수정하지 않았다. 단위격자 axis 부호/setting·기존 indexing의 확인이 필요하다.

조건부 B의 예측→metadata 등방 scale은 [10-1] DP01 약 0.950, [110] DM3 약 1.028, [001] TIFF 약 0.902, [110] TIFF 약 0.919이다. 서로 다른 배율/영상/입자에서 하나의 절대 calibration으로 정합되지 않는다. 시작 metric의 [001] fit은 배율 하한 0.900에 도달한다. 이 범위 내 최적값을 실측 calibration correction으로 적용하지 않았다.

추가 진단은 incremental 2×2 행렬 각 성분 ±0.04와 최종 중심 ±3 pixel을 허용한다. 실제 singular stretch는 별도 기록했다. 기준 parent에서 얻은 보정을 같은 그룹 모든 후보에 공유한다. 그 결과 DP01의 [10-1] RMS는 0.126→0.073, [001] TIFF는 0.085→0.041 nm⁻¹로 감소했다. [110] DM3는 0.147→0.142로 조금만 바뀌었고, [110] TIFF는 0.103→0.113으로 오히려 커졌다. Robust correspondence 목적함수와 trimmed RMS는 다르므로 모든 값의 단조 개선을 주장하지 않는다.

방향별 보정은 실제 격자 차이·검출기 왜곡·tilt·오지수화를 흡수할 수 있다. 이를 정련 격자나 독립 validation으로 취급하지 않는다. 원본과 CIF는 변경하지 않았다.

XRD 최상 `conditional_B_11b_seed29`의 기본 F1은 [10-1] DP01 0.141, [001] TIFF 0.452, [110] DM3 0.205, [110] TIFF 0.282다. 방향별 보정 진단에서는 각각 0.182/0.476/0.239/0.295다. 예를 들어 [001]의 ordered control은 같은 진단에서 0.750으로 더 높다. 모든 경우를 동시에 설명하는 구조로 채택할 근거는 없다. 이 수치를 정답 확률이나 실제 동역학 intensity fit으로 읽어서는 안 된다.

## q와 [110]에서 얻은 구별 조건

기존 DP01의 선택된 satellite-like 좌표 22개를 재사용한 상대 간격 비교다. 같은 영상의 재분석이므로 독립 증거를 추가한 것이 아니다.

| 가설 | 분수 간격 | ED 잔차 RMS / empirical v |
|---|---:|---:|
| 4b, m=1 | 0.250000 | 0.04018 |
| 7b, m=2 | 0.285714 | 0.00758 |
| 11b, m=3 | 0.272727 | 0.01812 |

ED 중앙값은 약 0.28823, XRD의 조건부 pure-b* 추정은 약 0.27826이다. ED 상대 간격은 2/7에 더 가깝지만 hkl·basis 및 공통 3D q는 미확정이다. 이 결과와 XRD의 11b 우세를 합쳐 주기를 임의 확정하지 않는다.

특히 **b 방향만 확장한 모든 후보의 정확한 [110] ZOLZ에서는 k_super=−n·h**이므로 parent k는 정수다. 4b/7b/11b 확장 자체로 이 평면에 fractional 위치를 추가할 수 없다. 점유 배열은 기존 위치의 intensity/소멸에는 영향을 준다. 실제 [110] 추가 반사가 별도 indexing·tilt/HOLZ·domain/중첩 검토 후에도 남는다면 c축 결합·적층·변위 등의 다른 가설이 필요하다. 이번 결과만으로 그 원인 하나를 확정하지 않는다.

## 산출물·검증

- `ed_matching_01/`: 모든 model/image/zone 반사 좌표·kinematic 강도·일대일 대응·threshold 민감도·회전 대조군·overlay.
- `ed_metric_diagnostic_01/`: 추가 방향별 보정 행렬 및 같은 후보의 보정 후 대응. 원본의 대체 결과가 아니다.
- `joint_review_01/`: 현재 process, q 비교, XRD/ED 점수 벡터와 한국어 보고서. 최종 shortlist는 비어 있다.

```bash
python -m lfnso.ed_matching configs/ed_matching.json
python -m lfnso.ed_metric_diagnostic configs/ed_metric_diagnostic.json
python -m lfnso.joint_review configs/joint_review.json
python -m pytest -q
```

원본 및 기존 추출 결과의 cloud 경로를 사용한다. 재실행은 새 output 경로로 수행한다. 39개 테스트가 통과했다. 새 검증에는 electron factor 교차 계산, 부분 점유 가중, exact repeat 소멸, zone 변환, 매칭 중복 방지·dense-grid precision, 합성 회전/배율 및 방향별 보정 복원을 포함한다. 이는 실제 구조 가설의 정답 검증이 아니다.
