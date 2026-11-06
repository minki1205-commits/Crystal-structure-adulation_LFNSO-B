# c축 가설 확장과 직접 XRD–ED 배열 탐색

2026-11-06 후속 작업. 최종 구조는 아직 미확정이다. 이번에는 c축 가설을 추가하고 **ED를 사후 비교뿐 아니라 탐색 목적함수에도 사용**했다. 다른 ED 방향에서 성능이 악화되는 현상까지 보존했다.

## 추가 계산

1. 기존 `starting_7b_seed17`, `conditional_B_7b_seed29`, `conditional_B_11b_seed29` 세 구조에서 출발했다.
2. Fe/Sb 교환 층의 primitive binary sequence를 2c/3c/4c/6c에서 순환 이동만 동치 처리해 열거했다. 각각 1/2/3/9개 sequence다. 역순 및 원소 교환을 동치로 간주하지 않았고, 실제 3D 구조의 완전한 symmetry-distinct 열거라고 주장하지 않는다.
3. Li/Sb/Fe-host 3상태 층별 교환, AABB 층 이동 b/48·b/24·b/6, cation breathing, 작은 층 이동과 Fe/Sb 교환의 결합 가설을 추가했다. Base 및 exact-repeat 대조를 포함해 **75개 시험 CIF**를 계산했다.
4. 이어 두 조건부 B base를 각각 6c로 확장하고 각 층 안에서 Li/Sb/Fe-host 자리를 직접 교환했다. XRD-only와 ED 가중치 0.4/1.2/3.0, seed 17/29를 비교해 **16개 초기 진단 종점 CIF**를 얻었다. 음수 background 허용의 영향을 확인한 뒤 비음수 background 조건으로 **16개를 다시 탐색**했다. 두 단계 합계 128,000 proposal / 86,744회 비동일 원소 교환 평가이며, 최종 비음수 단계는 64,000 proposal / 43,393회다. 반복 방문이 포함되므로 비등가 구조 개수와 같지 않다.

모든 모델은 정확한 Li4Fe0.67Ni0.33SbO6 평균 조성과 Fe/Ni .67/.33 공동 점유를 보존한다. 직접 탐색은 각 층의 Li:Sb:Fe-host 개수도 보존한다. 직접 탐색에서 산소/interlayer Li 좌표와 lattice metric은 고정했다. 에너지 relaxation, 산화수, 산소 결손은 적용하지 않았다.

Rigid registry 가설은 slab 산소와 slab 양이온을 함께 옮기고 interlayer Li에는 인접 층 이동의 평균을 적용했다. b/6 이동의 세 후보는 양이온 CN6 비율이 0.75여서 비교 적격에서 제외했다. 나머지 staged 72개와 직접 탐색 종점들은 해당 거리·점유·CN6 검사를 만족했다. **CN6는 cutoff 아래 이웃 개수 검사이며, octahedral bond angle·bond valence·에너지 안정성 증명이 아니다.**

Breathing은 slab Fe/Sb에 지정 진폭 0.04/0.08 Å, slab Li에는 −2배 진폭을 적용한 별도 변위 가설이다. 산소는 고정했고, 진폭은 측정값으로 제시한 것이 아니다.

## 직접 탐색의 정의

- XRD: 최종 단계는 20.4–34.6°에서 구조 scale과 선형 background 양끝을 비음수로 제한하고 R²의 제곱근을 최적화한다. 초기 단계만 unconstrained linear-background의 signed correlation을 사용했다. 잠정 Cu Kα 1.5406 Å, FWHM 0.25°.
- ED: [110] DM3와 기존 affine-diagnostic pose를 사용한다. 먼저 가능한 반사 격자와 관측 점의 일대일 대응을 정하고, simulated intensity를 최대값의 1%에서 부드럽게 포화시키는 **위치 support 점수**를 최적화한다. 실험 ED 강도 자체를 피팅하지 않는다.
- 최종 목적함수: sqrt(max(0, 비음수 background R²)) + 지정 가중치 × ED soft support. History의 기존 `XRD_correlation` 필드는 최종 단계에서 이 XRD 목적함수 성분을 기록한다. 가중치는 진단 설정이며 likelihood 또는 물리적 상대 오차에서 도출한 값이 아니다.
- 모든 최종 후보는 별도로 threshold 0.1/1/5%의 hard one-to-one matching, 기본 similarity pose와 affine diagnostic, 네 영상 전체에서 다시 평가했다.
- 다른 ED 방향 및 XRD 고각 구간은 목적함수에 직접 넣지 않았지만 이전 탐색에서 이미 본 자료이므로 새로운 독립 검증으로 부르지 않는다.
- 부모 구조의 배율·회전·왜곡 진단을 고정해 공유했다. 후보에 맞춰 ED calibration을 다시 조정하지 않았다.

## 배경 제약 재검토

초기 최고 R² 0.849는 음수 선형 background를 허용한 결과다. 같은 CIF를 비음수 background로 다시 평가하면 약 0.783으로 내려간다. 이 값을 최종 개선 성과로 사용하지 않고 제약을 탐색 자체에 넣어 재계산했다. **이후 새 원자 배열에서 얻은 최고 값은 0.846**이며 아래 표는 이 비음수 재탐색 기준이다. 최고 후보는 background 양끝이 0에 도달하므로, 기기 배경·peak width·많은 약한 반사의 기여를 분리하지 못한 점도 남아 있다. 이 제약만으로 물리적 background가 검증된 것은 아니다.

## 무엇이 개선되고 무엇이 남았는가

아래 ED는 affine diagnostic, threshold 1%, DM3 위치 허용오차 0.08 nm⁻¹ 조건이다. DP01은 미확정 [10-1] 대안 setting으로 평가했다.

| 설명용 후보 | XRD superstructure R² | [110] DM3 F1 | DP01 F1 |
|---|---:|---:|---:|
| 이전 조건부 B 11b | 0.705 | 0.239 | 0.182 |
| 11b×6c, XRD-only, seed29 | 0.846 | 0.239 | 0.099 |
| 11b×6c, ED weight 3.0, seed17 | 0.737 | 0.331 | 0.099 |
| 7b×6c, ED weight 3.0, seed29 | 0.612 | 0.360 | 0.115 |

정해진 층 순서·변위만 시험했을 때의 최고 XRD 값은 약 0.711로 개선이 작았다. 비음수 배경 조건으로 층별 배열을 직접 바꾼 뒤 11b×6c에서 0.846까지 올라갔다. 같은 후보의 35–45° R²는 약 0.789이다. 그러나 [101]-named 영상의 불일치가 커졌다.

ED 비중을 늘리면 [110] DM3/TIFF 대응이 개선되지만 XRD 점수가 내려간다. 7b×6c의 ED-weight=3 후보는 [110] TIFF F1도 약 0.444로 증가했다. 그래도 [001]/DP01을 포함한 공동 검증을 통과하지 않았다. 한 방향의 성능이나 두 목적함수의 합으로 최종 승자를 선언하지 않았다.

최종 R²는 비음수 linear-background-only fit 대비 SSE 개선량으로 Rwp나 정답 확률이 아니다. F1은 위치 대응의 precision/recall 조화평균으로 동역학적 전자 강도 적합도가 아니다. 6c는 이번에 시험한 근사 주기이며 발견된 실제 주기로 보고하지 않는다.

## [110] 위치 격자의 상한 진단

추가 구조 전에 intensity를 무시한 c-repeat 격자 1/2/3/4/6/7/8/12를 비교했다. DM3의 87개 점에 대해 기본 c격자는 14개, 6c 격자는 43개를 허용오차 내에서 대응할 수 있었다. 그러나 예측 위치 수도 91→557개로 늘었다. 12c에서는 51개를 맞추면서 1,113개를 예측한다. 따라서 grid coverage만으로 6c/12c를 선택할 수 없다. 모든 숫자는 기존 pose·선택 ROI·threshold 아래의 조건부 진단이다.

## 파일과 검증

- `coupled_search_01/`: 75개 CIF, metadata, XRD 및 네 영상의 반사표, 모든 매칭 조건.
- `joint_occupational_01/`: ED weight 0/0.4 직접 탐색 8개.
- `joint_occupational_02/`: ED weight 1.2/3.0 직접 탐색 8개.
- `joint_occupational_03/`: 비음수 background 조건으로 다시 탐색한 최종 비교용 16개.
- `coupled_review_01/`: 초기 unconstrained-background 진단 기록. 최종 비교에 사용하지 않는다.
- `coupled_review_02/`: 최신 보고서, 비음수 XRD overlay, tradeoff 그림, score CSV, 기존 음수배경 결과의 재평가, 대표 CIF의 독립 검증.
- `results/LFNSO_3D_joint_search_results.zip`: 검토용 결과 묶음. 원본 현미경 파일은 포함하지 않는다.
- `results/latest_joint_review.html`: 이번 최신 보고서 별칭. 이전 보고서는 앞선 ED 결과 ZIP과 원래 cycle에 보존했다.

대표 세 CIF를 다시 읽고 독립 pymatgen powder 계산으로 검증했다. 상대 강도 최대 차이는 약 3.6×10⁻¹² 이하였다. 이는 구현·파일 입출력의 수치 검증이며 실험 구조의 정답 검증이 아니다. Exact-repeat 대조와 layer stoichiometry 검사, 3D zone law, 산소·interlayer 처리, cache swap 및 checkpoint 검사 등을 포함한 **52개 테스트가 통과**했다.

```bash
python -m lfnso.coupled_search configs/coupled_search.json
python -m lfnso.joint_occupational configs/joint_occupational.json
python -m lfnso.joint_occupational configs/joint_occupational_strong_ed.json
python -m lfnso.joint_occupational configs/joint_occupational_nonnegative.json
python -m lfnso.coupled_review configs/coupled_review_nonnegative.json
```

새 실행에서는 output 경로를 새 cycle로 지정한다. 직접 탐색은 동일 입력·코드 signature의 checkpoint에서 재개할 수 있다. 각 cycle의 실행 당시 source snapshot을 보존했다. 원본 경로 및 이전 분석 의존성은 cloud 설정 파일에 명시되어 있다.

## 남은 핵심 제약

현재 가장 큰 문제는 c축 자유도로 [110]/XRD를 개선하면서 in-plane satellite ordering의 ED 정합을 유지하지 못한다는 점이다. 다음 탐색은 [001]·[101]/[10-1] 위성 조건도 함께 보존하는 다중-zone 제약, 검증된 축 convention/calibration, 제한된 metric/oxygen 변위의 영향을 분리해야 한다. 현재 결과로 최종 구조, 특정 적층 순서 또는 STEM 검증 통과를 확정하지 않는다.
