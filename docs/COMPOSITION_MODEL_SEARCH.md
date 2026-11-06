# Li4Fe0.67Ni0.33SbO6 가설 구조 생성

사용자가 제공한 조성 및 Ni의 Fe-site 치환 가정을 반영하여, 기존 Li2MnO3 CIF 좌표에서 생성 작업을 진행했다. 기존 LFNSO CIF가 없어도 구조 가설을 만들 수 있도록 site mapping을 명시했다. 이 mapping은 검증된 LFNSO 원자 배열이 아니다.

## Parent와 배위

제공된 conventional Li2MnO3 셀은 Li8Mn4O12이다. Li/O 좌표는 보존하고 네 Mn site 중 두 개에 Sb, 두 개에 Fe0.67Ni0.33를 배정하면 Li8Fe1.34Ni0.66Sb2O12, 즉 사용자 조성의 두 formula unit이 된다. Ni와 Fe는 같은 host의 점유 비율을 유지한다. Li의 slab/interlayer 구분은 실제 fractional z 좌표에서 얻었다.

원본에서 Li와 Mn의 2.6 Å 이내 O coordination은 모두 6이었다. 시작 LFNSO metric으로 바꾼 parent와 대다수 후보에서도 2.7 Å cutoff의 sixfold coordination을 확인했다. 산소 좌표는 Li2MnO3에서 상속한 값으로 에너지 relaxation은 하지 않았다.

Parent symmetry로 구별되는 Sb 배정은 3개이며 원본 site index 쌍은 (8,9), (8,10), (8,11)이다. P1 CIF 표현은 계산용이며 실제 공간군 판정이 아니다. Parent 격자는 master prompt의 a=5.1706, b=8.9382, c=5.1635 Å, β=109.49°로 고정했다.

## 생성한 32개

- 서로 다른 Sb parent 배열 3개.
- 3b, 4b, 6b, 7b, 11b, 12b에서 exact-repeat, Li/Fe-Ni 교환, coherent displacement 모델 18개.
- AABB stacking 1개와 occupancy+stacking coupled 모델 2개.
- Fe/Ni 정수 원자 배열 근사 1개.
- 7b의 Sb/Fe-Ni host 교환 3개, 더 큰 Li/TM 교환·변위 진폭 4개.

7b는 2/7 b*, 11b는 3/11 b*를 시험하는 rational approximant이다. 물리적으로 정확한 주기라는 주장이나 관측 q를 임의 확정한 것이 아니다. 12b의 m=3 변조는 4b 반복을 갖는 대조군이다. 단순 supercell 확장과 실제 변조를 구분한다.

31개 평균 점유 모델은 정확히 Li4Fe0.67Ni0.33SbO6를 보존한다. 정수 배열 1개는 Fe2/3Ni1/3로 반올림한 **명시적 근사**이며 정확한 .67/.33 모델과 동일 취급하지 않는다. 부분 점유 CIF는 평균 산란 구조이고 하나의 구체적인 원자 배열이 아니다.

AABB는 각 TM slab의 O 및 slab Li까지 함께 이동시켰고 interlayer Li는 인접 slab 이동의 평균을 적용했다. 그런데 이 시험 배치에서는 cation의 25%가 2.7 Å cutoff 아래 sixfold가 아니었다. 해당 세 모델은 구조 검토 대상으로 남기고 검증 후보로 채택하지 않았다.

## 회절과 실제 데이터 비교

- 32개 모두 17–36°의 ideal periodic XRD를 계산했다. Cu Kα 1.5406 Å는 기존 잠정 가정이다.
- [001]/[101]/[110] 기하학적 반사표와 Z-weight structure-factor proxy를 출력했다. **Z proxy는 실제 electron scattering factor나 dynamical ED intensity가 아니다.**
- Exact-repeat 6개에서 parent와 XRD 위치 및 상대 강도가 동일한 것을 확인했다.
- 실험의 다섯 feature에 일대일 peak matching을 적용했다. 검출 threshold .1/.5/1%는 분석 민감도 설정이며 실제 experimental detection limit이 아니다.
- 20.4–34.6° 실험 수치 프로파일과도 비교했다. Gaussian width .15/.25/.4°를 시험하고, 각 모델에 하나의 양수 intensity scale과 선형 background만 피팅했다. 이 trial은 Rietveld refinement가 아니다.

7b의 Sb/Fe-host pattern1 및 최대 Li/TM 교환 모델은 선택된 feature 위치의 RMS를 약 0.0576°까지 줄였다. 그러나 이 모델들도 전체 해당 구간에서 B 강도 분포를 설명하지 못했다. 현재 32개 모두 profile fit의 nonnegative structural scale이 사실상 0으로 내려가고 background-only 수준에 머물렀다. 즉 peak 위치 일부를 맞추는 것과 실제 패턴을 설명하는 것은 다르다.

따라서 **검증된 최종 shortlist는 없다**. 수치 반올림 차이로 exact-repeat 모델을 서로 다른 순위로 취급하지 않도록 Pareto 비교에 tolerance를 적용했다. 강도 설명력이 없는 모델을 Pareto 채택으로 오인하지 않도록 별도 flag를 기록했다.

## 조성·전하와 남은 연구

조성과 배위만으로 산화수는 정해지지 않는다. Li(+1), Fe(+3), Sb(+5), O(-2)를 가정할 때 Ni(+2) 시 formula당 -0.33e, Ni(+3) 시 0이다. 이는 검토할 시나리오일 뿐이며 Ni 산화수나 산소 vacancy를 임의로 부여하지 않았다.

현재 결과는 이 제한된 배열·진폭·고정 metric 가설들의 실패를 보여준다. 모든 occupational/stacking/displacement 가설이 배제된 것은 아니다. 다음 탐색은 Fe/Ni/Sb/Li의 더 다양한 slab 내 배열과 평균 metric/oxygen relaxation을 구분하고, 강한 불필요 반사를 억제하면서 B의 주요 반사를 만드는 배열을 찾는 데 초점을 둬야 한다. XRD 및 다중 zone gate를 통과한 구조가 없어 STEM multislice는 실행하지 않았다.

## 파일 및 재현

- 기본 25개: `03_structure_model_search/composition_models_02/`
- 진폭·host-exchange 확장 7개: `03_structure_model_search/composition_models_03/`
- 각 후보: `model.cif`, `metadata.json`, `xrd.json`, `zone_*.csv`
- 각 cycle: `summary.json`, `assessment.json`, `screening.csv`, `pareto_scores.csv`, `report_assessed.html`

실행 명령은 저장소 루트에서 `python -m lfnso.model_search configs/model_search.json` 및 `configs/model_search_extension.json`이다. Python은 `/workspace/lfnso-venv/bin/python`, 환경은 `PYTHONPATH=src`, `MPLCONFIGDIR=/workspace/lfnso-cache/matplotlib`, `OPENBLAS_NUM_THREADS=1`을 사용한다. 새 cycle output을 지정해야 하며 기존 결과를 덮어쓰지 않는다.

초기 개발 중 중단된 `composition_models_01`과 `assessment_initial_numeric_tie_not_for_interpretation.json`은 최종 결과가 아니다. 배포 ZIP에서 제외했다. 원본 CIF와 실험 데이터는 수정하지 않았다.
