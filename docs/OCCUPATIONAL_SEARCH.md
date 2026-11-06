# 조성 보존 자리 배열 탐색 — 2026-11-05

Li4Fe0.67Ni0.33SbO6 조성, Ni의 Fe-host 치환, 팔면체 배위라는 사용자 조건을 반영해 **20개 탐색 종점 CIF**를 추가했다. 산화수·산소 결손은 부여하지 않았다. 기존 32개 후보 및 원본 자료는 보존했다. 이번 결과도 검증 전 가설이며 최종 구조 shortlist는 없다.

## 무엇을 바꾸었는가

기존 Li2MnO3 좌표의 slab 내 여섯 양이온 자리에 Li, Sb, Fe0.67Ni0.33를 교환 배치했다. conventional parent당 각 종류 두 자리의 총량을 supercell 전체에서 보존한다. 각 tile의 국소 조성은 달라질 수 있다. Interlayer Li와 O 좌표는 고정했다. Fe/Ni는 항상 같은 host의 .67/.33 평균 점유이며 개별 Ni 원자 배열은 정하지 않았다.

4b/7b/11b 셀은 각각 96/168/264 site이다. 모든 종점에서 정확한 평균 조성, 자리 점유 합 1, 2.7 Å cutoff 아래 양이온의 O 이웃 6개를 확인했다. 이는 물리적 안정성이나 결합 원자가 검증이 아니다. 산소 좌표는 상속한 값이며 에너지 relaxation을 수행하지 않았다.

두 metric을 비교했다.

- `starting`: a=5.1706, b=8.9382, c=5.1635 Å, β=109.49°.
- `conditional_B`: a=5.1543501, b=8.9454194, c=5.1437185 Å, β=108.8042271°. 앞선 다섯 피크·다섯 변수 조건부 해이며 자유도 0으로, 독립 정련 결과가 아니다.

## 계산 및 재현성

각 시작 배열에서 조성을 보존하는 두 자리 교환을 수행하고 simulated annealing으로 프로파일과의 signed correlation을 최적화했다. 후보 XRD는 Cu Kα **잠정 파장 1.5406 Å**, 고정 Gaussian FWHM 0.25°로 계산했다. 0.15/0.40°는 같은 배열의 사후 민감도 검사이며 그 폭에서 재탐색한 결과가 아니다.

- 1차: 두 metric × 7b/11b × seed 17/29, 각 2,000 proposal. 20.4–30.2°만 최적화하고 30.3–34.6°와 35–45°를 사후 비교.
- 2차: 두 metric × 4b/7b/11b × seed 17/29, 각 6,000 proposal. 34.35° 실패를 확인한 후 20.4–34.6°까지 탐색 구간을 확대.
- 합계 **88,000 proposal**, 동일 원소 쌍을 제외한 **60,230회 회절 목적함수 평가**. 반복 방문을 포함하므로 서로 다른 구조 60,230개라는 뜻이 아니다. 대칭 비등가 구조의 완전 열거 또는 전역 최적해를 주장하지 않는다.
- Checkpoint에는 현재/최상 배열, RNG 상태, 진행 단계, 입력·코드 signature를 저장한다. 중단 후 이어서 실행할 수 있고 다른 signature는 거절한다.
- 교환으로 갱신하는 structure factor를 pymatgen XRD 계산과 교차 확인했다. 모든 최종 CIF를 다시 읽어 회절을 확인했다. 재독출에는 `CifParser(..., frac_tolerance=0)`을 사용해 reader의 자동 유리수 좌표 변환을 끈다.
- 테스트는 구조 인자 교차 확인, 100회 교환 후 전체 재계산 일치, 주기 경계 slab 분류, 조성 보존, 합성 프로파일 scale/background 복원, 실제 중단을 모사한 checkpoint 재개, Fourier 주기 구분을 포함한다.

## 결과

아래 R²는 `1 − SSE(model+linear background)/SSE(linear background only)`이다. **Rwp·구조 정답 확률·설명된 원자 비율이 아니다.** 음수 구조 scale을 허용하지 않으며, 각 구간의 scale/background를 따로 맞췄다. 고각 점수는 절대 강도 전달 검증이 아니다.

| 2차 탐색 조건 | 20.4–34.6° R², 두 seed 범위 | 35–45° R², 두 seed 범위 |
|---|---:|---:|
| 시작 격자 · 4b | 0.275–0.278 | 0.489–0.495 |
| 시작 격자 · 7b | 0.654–0.656 | 0.502–0.502 |
| 시작 격자 · 11b | 0.582–0.588 | 0.503–0.506 |
| 조건부 B 격자 · 4b | 0.223–0.225 | 0.773–0.777 |
| 조건부 B 격자 · 7b | 0.576–0.583 | 0.762–0.767 |
| 조건부 B 격자 · 11b | 0.696–0.705 | 0.753–0.761 |

1차 탐색은 선택 구간에서 R² 0.646–0.768을 얻었지만, 남겨 둔 30.3–34.6°에서는 0–0.014에 그쳤다. 2차의 `conditional_B_11b_seed29`는 20.4–34.6° 점수가 가장 높고 34.35° 부근에도 반사를 만든다. 그러나 24.8° 부근의 모양·강도 불일치와 30–33°의 불필요한 반사들이 남는다. 30.3–34.6°만 별도 평가하면 R²는 0.225에 불과하다.

조건부 B 격자의 ordered control은 35–45°에서 R²≈0.773이다. 최고 superstructure 후보의 같은 구간 값 0.753보다 높다. 따라서 고각 개선을 새 ordering의 증거로 해석할 수 없고, 격자 효과를 구분해야 한다. 후속 설계 전에 이미 이 구간을 확인했으므로 2차 결과를 새 독립 검증으로 부르지 않는다.

피크 폭 민감도도 크다. 최고 후보의 20.4–34.6° R²는 FWHM 0.15/0.25/0.40°에서 0.569/0.705/0.667, 저각 17–20.3°에서는 0.910/0.783/0.577이다. 하나의 고정 폭만으로 전체 데이터를 평가하는 한계가 보이며, 이를 구조 확정으로 보완해서는 안 된다.

7b 및 11b는 여러 harmonic을 허용한다. 이번 2차 후보의 slot별 occupational Fourier power는 7b에서 대체로 m=2, 11b에서 m=3가 가장 크지만 다른 성분도 상당하다. 최고 11b 후보의 m=3 성분은 non-DC power의 약 35.4%이다. 이는 순수 q=3/11 b* 변조나 ED 위성 강도가 확인됐다는 뜻이 아니다. Fourier power는 원소 종류를 one-hot 표시한 배열의 기술 통계이다.

## 판단 및 다음 계산의 기준

이번 탐색에서는 slab 자리 배열 확장으로 기존 모델의 강도 실패를 일부 개선했다. 단순한 4b 가설보다 7b/11b 후보가 이 제한된 superstructure 목적함수에서 우세하지만, 주기·원자 배열은 미확정이다. 모든 후보의 XRD/다중-zone gate 통과가 성립하지 않아 STEM multislice를 실행하지 않았다.

다음 단계에는 먼저 격자·zero offset·피크 폭/비대칭과 ordering의 영향을 분리한 전체 XRD 비교가 필요하다. 산소 및 양이온 변위·적층 가설은 각각 독립된 제한 모델로 추가해야 한다. ED의 축 convention/calibration과 다중-zone 위치 제약도 해결해야 한다. 이 단계에서 사용자의 추가 자료 없이 산화수, 산소 결손, 두 상 공존, 정확한 q를 확정하지 않는다.

## 산출물 및 실행

- `03_structure_model_search/occupational_search_02/`: 1차 8개 종점, CIF, metadata, XRD, 비교 그림, checkpoint.
- `03_structure_model_search/occupational_search_03/`: 2차 12개 종점과 동일 산출물.
- `03_structure_model_search/occupational_assessment_03/`: 최종 진단 보고서, 전체 점수 CSV/JSON, 그림, 실행 소스 snapshot.
- `03_structure_model_search/LFNSO_occupational_search_results.zip`: 위 결과와 재현 코드/설정/설명.
- `occupational_search_01`, `occupational_assessment_01`은 개발 중 중단된 디렉터리이며 결과로 사용하지 않는다. `occupational_assessment_02`는 앞선 동일 수치 진단이며 최신 보고서는 03이다.

저장소 루트에서 준비된 venv로 실행한다. 완료된 탐색은 checkpoint에서 재개하며 추가 proposal 없이 검증만 수행할 수 있다. 새 입력·코드·탐색 조건을 사용할 때에는 새 output 경로를 지정한다. Assessment는 기존 output 경로에 덮어쓰기를 거절한다.

```bash
python -m lfnso.occupational_search configs/occupational_search.json
python -m lfnso.occupational_search configs/occupational_search_full_superstructure.json
python -m lfnso.occupational_assessment configs/occupational_assessment.json
python -m pytest -q
```

원본은 `/workspace/lfnso-data/`에 있고 설정에는 현재 cloud 경로가 기록되어 있다. 다른 환경에서는 원본 경로를 조정해야 한다. ZIP은 이번 계산의 검토용 묶음이며 전체 microscopy 원본을 포함하지 않는다. Cloud workspace 저장과 원격 Git 동기화는 별개이며 이 단계에서 원격 push를 수행하지 않았다.
