# LFNSO21-B structure discovery

Li4Fe0.67Ni0.33SbO6의 XRD·ED·STEM 자료를 함께 설명할 구조 가설을 탐색하는 연구 코드입니다. 자체 생성 모델은 Ni를 Fe-host에 .67/.33 평균 점유로 배치하고, 공동연구기관 CIF는 제공된 층별 혼합 점유 그대로 별도 검토합니다. 현재 **후보 생성 → XRD screening → 실제 ED 위치 정합**까지 진행했고, 공동 검증을 통과한 최종 구조는 없습니다. STEM multislice는 아직 실행하지 않았습니다.

## 최신 결과

- [공동연구기관 C2 P/S 교대 CIF 검토](docs/COLLABORATOR_CIF_REVIEW.md): 받은 CIF는 XRD 전체를 설명하지 못하지만 [110] ED와 [001] 기본 격자에는 부분 대응이 있습니다. 저각 피크에 격자를 맞추면 35–45° 구간이 악화되며 DP01 대응도 미해결입니다. 사용자가 알려준 이번 CIF의 STEM [001]/[100]과 ED 파일 축을 구분해 비교했습니다. 보고서: `results/collaborator_CIF_review.html`.
- [c축 가설 확장 및 XRD–ED 직접 배열 탐색](docs/COUPLED_SEARCH.md): 75개 층별 시험 구조와 32개 직접 탐색 종점(초기 진단 16개·비음수 배경 재탐색 16개)을 추가했습니다. 비음수 배경 조건의 XRD-only 최고 R²는 0.846으로 개선됐지만 다른 ED 방향의 불일치가 남아 최종 구조는 미확정입니다. 최신 보고서는 `results/latest_joint_review.html`입니다.
- [현재 process 및 실제 ED 비교](docs/ED_MATCHING.md): 후보 12개 + 대조 구조 2개, DM3/TIFF 네 영상, zone 대안 포함 70개 조합. 위치·반사 threshold·배율/왜곡 민감도와 미설명 반사를 기록했습니다.
- [조성 보존 배열 탐색](docs/OCCUPATIONAL_SEARCH.md): 4b/7b/11b와 두 metric, 탐색 종점 CIF 20개. XRD 일부 구간의 11b 우세는 최종 주기 확정을 뜻하지 않습니다.
- [이전 32개 가설](docs/COMPOSITION_MODEL_SEARCH.md), [수치 XRD 분석](docs/RAW_XRD_ANALYSIS.md), [초기 위치 가설](docs/POSITION_SCREENING.md), [입력 자료 점검](docs/FIRST_DATA_SURVEY.md).
- [원래 연구 계획](docs/MASTER_PROMPT.md). 초기 분석 문서는 당시 상태를 보존한 기록이므로 현재 완료 범위는 최신 문서를 기준으로 봅니다.

`results/`에는 검토용 결과 묶음과 현재 공동 검토 보고서를 보관합니다. 새 계산의 전체 산출물은 Git에서 제외된 `03_structure_model_search/`에 생성합니다. 원본 현미경 파일은 별도 보관하며 현재 cloud 경로는 `/workspace/lfnso-data/`입니다. 검토용 묶음은 원본 현미경 데이터 전체의 백업이 아닙니다.

## 실행

Python 3.12, 저장소 루트 기준:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock.txt
.venv/bin/python -m pip install --no-deps -e .
.venv/bin/python -m pytest -q
.venv/bin/python -m lfnso.ed_matching configs/ed_matching.json
.venv/bin/python -m lfnso.ed_metric_diagnostic configs/ed_metric_diagnostic.json
.venv/bin/python -m lfnso.joint_review configs/joint_review.json
```

Cloud 전용 설치는 `bash scripts/setup.sh`를 사용합니다. 분석 config는 현재 cloud 원본/이전 결과 경로를 참조하므로 다른 환경에서 원본 경로를 맞춰야 합니다. 각 output은 새 cycle을 지정하며 기존/원본을 덮어쓰지 않습니다. Occupancy 탐색은 일치하는 입력·코드 signature의 checkpoint에서 재개할 수 있습니다. 검토 ZIP의 `REPRODUCE.txt`에는 해당 묶음의 별도 재현 절차가 있습니다.

## 해석의 범위

Cu Kα 1.5406 Å는 사용자 제안에 따른 잠정 가정입니다. 후보는 Li2MnO3에서 상속한 산소 좌표를 사용한 미정련 구조이며, 조성·배위 검사는 에너지 안정성 검증이 아닙니다. 산화수나 산소 결손을 확정하지 않았습니다.

ED는 전자 산란인자를 사용한 단일 산란 근사와 위치 매칭입니다. 실험 강도 피팅·두께/tilt/다중 산란·최종 공간군 판정은 하지 않았습니다. [101]/[10-1] setting과 영상 calibration이 미확정이며, 방향별 보정 결과는 민감도 진단입니다. 저장 FFT는 원본 HAADF와 중복된 증거로 취급합니다.

## 직장·집에서 작업 이어가기

작업 브랜치는 `Structure-model-adulation`입니다. 기존 checkout에서 해당 브랜치를 pull한 후 같은 환경 설정과 원본 자료 경로로 작업합니다. Git 코드 동기화와 cloud 환경 설정 게시, 원본 데이터 백업은 각각 별개입니다. 새로운 결과는 코드·설정·출처 hash와 함께 검토 가능한 묶음으로 보존합니다.
