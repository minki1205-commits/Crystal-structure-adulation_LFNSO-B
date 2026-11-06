# 첨부 데이터 1차 분석

확인된 정보: 사용자에 따르면 모든 첨부 영상은 LFNSO21-B이며 image1은 HAADF이다. Reference/A 영상은 없다. 파일명 zone axis는 사용자 표기를 보존했으며 독립 지수화는 아직 수행하지 않았다.

## 확보 및 원본 보존

- 원본 복사본: `/workspace/lfnso-data/raw`, 출처·SHA256: `/workspace/lfnso-data/manifest.json`.
- 11개 첨부 중 CIF 1개와 DM 9개 확보. `[001] OneView_300X_200kV_0021_001.dm4`는 전송 도구의 32 MiB 제한으로 미확보.
- Windows 기존 작업 폴더는 클라우드에서 접근할 수 없다.
- 기존/원본 파일을 수정하지 않았다. 출력은 새 cycle 디렉터리에 저장한다.

## 관측·검증 결과

- CIF는 C2/m Li₂MnO₃, conventional cell Li8Mn4O12, 24 sites. a=4.937 Å, b=8.532 Å, c=5.030 Å, β=109.46°. LFNSO의 화학적 parent가 아닌 기하학적 reference이다. Li(+1), Mn(+4), O(-2)를 가정하면 중성이나 이는 측정 산화수 검증이 아니다.
- 이 CIF의 1b/3b/4b/6b 단순 반복 대조 구조 4개와 회절 좌표표를 생성했다. 실제 LFNSO 후보로 순위를 매기지 않았다. 파장이 제공되지 않아 이번 실제 입력에 대한 XRD는 실행하지 않았다.
- DP01_101 및 dp3_110은 2048², 0.016421644 nm⁻¹/pixel. 69개/107개 자동 피크 후보를 추출했다. 알고리즘 threshold에 의존하며 artifact·중복·누락 가능성이 있다.
- 반대편 반사쌍을 통한 중심 추정: [101] y=820.695, x=935.215 px (23쌍), [110] y=774.199, x=841.880 px (21쌍). 이 과정은 intensity Friedel symmetry를 가정하지 않는다. 좌표축은 영상 y-down/x-right이며 아직 결정축에 회전 정합하지 않았다.
- 20MX 영상: 0.009347768 nm/pixel, 시야 약 9.572 nm, FFT bin 약 0.104470 nm⁻¹.
- 12MX 영상: 0.015579613 nm/pixel, 시야 약 15.954 nm, FFT bin 약 0.062682 nm⁻¹.
- LFNSO 프롬프트의 b를 사용한 가설 q3≈0.37293, q4≈0.27970 nm⁻¹의 차이는 약 0.09323 nm⁻¹이다. 20MX 한 영상의 nominal Fourier bin보다 작다. Subpixel fitting을 원천적으로 배제하는 뜻은 아니지만, window leakage·noise·scan distortion을 무시한 단순 bin 판독으로 두 가설을 구별할 수 없다.
- 저장 FFT와 HAADF 0075의 재계산 unwindowed FFT는 DC 제외 진폭 상대 L2 오차 약 1.43×10⁻⁶로 일치한다. 복소 위상 convention은 다르므로 위상 재구성에는 원본 영상 FFT를 사용해야 한다. 독립 증거로 중복 집계하지 않는다.
- BF/HAADF 0076 및 0077은 촬영 시각이 일치하는 paired channels로 기록했다. 두 촬영 간 같은 입자/영역 여부는 확인되지 않았다.
- OneView_250X_200kV_0052_110.dm4를 읽은 배열에는 NaN 16개와 float 범위에 가까운 값이 존재한다. 원본 손상인지 export/reader 문제인지 확정하지 않았으며 정량 분석에서 격리했다. 임의 보정하지 않았다.

## 재현

저장소 루트에서 `bash scripts/setup.sh` 후 다음을 실행한다. 기존 cycle이 있으면 config의 output을 새로운 이름으로 바꾼다.

```bash
export MPLCONFIGDIR=/workspace/lfnso-cache/matplotlib
PYTHONPATH=src /workspace/lfnso-venv/bin/python -m lfnso.survey configs/uploaded_survey.json
PYTHONPATH=src /workspace/lfnso-venv/bin/python -m lfnso.pipeline configs/li2mno3_controls.json
```

결과: `03_structure_model_search/uploaded_survey_02/` (HTML/JSON/CSV/미리보기), `03_structure_model_search/li2mno3_controls_01/` (대조 CIF/회절 좌표).

## 다음 입력의 우선순위

1. 기존 LFNSO reference CIF와 4b/6b/AABB 생성 코드·report. Li/Fe/Ni/Sb site assignment와 조성 확인에 필요하다. Li₂MnO₃ 좌표에 임의 치환을 적용하지 않는다.
2. XRD raw profile과 파장, 시료 식별, 측정 조건. Gate 1과 평균 격자 변화/ordering 반사 비교에 필요하다.
3. 미확보 [001] OneView의 32 MiB 이하 분할 압축 또는 TIFF+calibration, 그리고 이상값이 없는 [110] OneView 재export. TIFF에는 pixel size, 단위, 원본 크기와 회전/crop 여부를 함께 제공한다.
4. 기존 [110] FFT 분석 report와 사용한 crop/mask, hkl 표기 자료. HAADF 수렴각·검출기 내/외각·두께 정보는 이후 STEM 단계에서 필요하다.

현재 XRD/ED gate는 not evaluated, structural ranking은 미실행이다. 피크 좌표 추출은 구조 판정이 아니다.
