# XRD 수치 프로파일을 이용한 갱신 분석

## 입력과 파장 가정

`LFSO_LFNSO21 XRD data.xlsx`의 `그래프데이터` 시트에서 LFSO / LFNSO21-A / LFNSO21-B를 읽었다. 각각 2θ=10–80°, step=0.02°, 3501 points이며 A:B, C:D, E:F 열이 angle/intensity 쌍이다. 결측·비수치·비단조 각도는 reader에서 거부한다.

원본 복사본은 `/workspace/lfnso-data/raw_xrd/`, SHA256과 업로드 출처는 `/workspace/lfnso-data/xrd_manifest.json`에 보존했다. 데이터에 새로운 평활화나 배경 제거를 적용해 저장하지 않았다. 전처리 이력은 엑셀에 없고 사용자가 확인하지 못했으므로 `workbook intensity as supplied`로 취급한다. 독립 Poisson detector counts라고 가정하지 않는다.

사용자가 MiniXRD이므로 Cu Kalpha로 가정하도록 제안했다. **λ=1.5406 Å는 잠정 기준**이며 확인된 기기 설정이 아니다. λ=1.5418 Å sensitivity도 기록했다. Kα₂ 제거/단색화, 표준시료 보정과 zero shift는 미확인이다.

## 피크 위치 측정

22개 주요 feature에 대해 평활화하지 않은 수치값에 pseudo-Voigt envelope + 국소 선형 배경을 피팅했다. Half-window 0.25/0.35/0.45°로 중심 민감도를 확인했다. 출력의 formal SE는 local least-squares noise 가정 아래 값이며 기기·Kα·전처리·overlap·shape model 오차를 포함하지 않는다. 관측 FWHM을 결정 크기/결함 밀도로 바로 변환하지 않는다.

| B feature | 2θ / ° | d / Å (잠정 Cu Kα) |
|---|---:|---:|
| (001)로 가정한 anchor | 18.2048 | 4.86918 |
| (020)로 가정한 anchor | 19.8341 | 4.47271 |
| 추가 feature 1 | 22.1924 | 4.00245 |
| 추가 feature 2 | 24.7602 | 3.59288 |
| 추가 feature 3 | 29.1943 | 3.05648 |
| 추가 feature 4 | 34.3535 | 2.60836 |

이 수치가 앞선 그림 수동 판독값을 대체한다. 특히 세 feature의 그림 판독 4.03/3.62/3.08 Å는 현재 잠정 파장 변환값보다 약 0.028/0.027/0.024 Å 크다. 따라서 그림 기반으로 평가한 1/4 b*의 높은 일치도를 원본 수치의 정밀 일치로 취급하면 안 된다. 파장 1.5418 Å를 택할 때 d는 약 0.078% 증가하며 이 차이를 모두 설명하지는 않는다.

A의 주요 feature는 20.7773°, 23.2523°, 27.0980°, 31.7230°이다. 모든 A/B feature가 동일 hkl이라는 가정을 먼저 검증해야 한다.

## 파장을 소거한 조건부 변조 비율

Unique-b monoclinic setting에서 B의 세 반사를 (001), (020), (0,2+δ,1)로 대응시키면, s=2 sinθ에 대해

`δ = 2 sqrt[(s_sat² − s_001²) / s_020²] − 2`

이다. 이 관계는 파장·a·β에 의존하지 않는다. 다만 **선택한 hkl 대응과 순수 b* 변조 가정**에 의존한다.

- δ≈0.27826; local formal propagated SE≈0.00039. 이 작은 formal SE를 최종 구조 불확실성으로 해석하지 않는다.
- 피팅 window 변화 결과는 약 0.27794–0.27860.
- 진단용 zero-shift ±0.05°를 적용하면 약 0.27484–0.28169. 이 ±0.05°는 확인된 기기 오차 범위가 아니라 민감도 시험이다.
- 같은 두 anchor에서 δ=1/4를 강제하면 세 번째 feature 예측은 28.9707°, 관측은 29.1943°로 약 0.224° 차이가 난다.
- δ=1/3의 예측은 29.6338°로 약 0.439° 차이가 난다.

그러므로 이 반사 대응 아래에서 **정확한 4b 또는 3b를 확정하기 어렵다**. 자유 q, 다른 indexing, 다른 성분의 q 또는 반사족 대응을 계속 검증해야 한다. 많은 rational approximant를 추가해 숫자에 맞추는 방식으로 주기를 확정하지 않는다.

## ED와의 연결

잠정 Cu Kα와 (020) 대응에서 b*=1.11789 nm⁻¹, q=δb*≈0.31106 nm⁻¹, 1/q≈3.215 nm이다. 이전 DP01_101의 경험적 q≈0.30940 nm⁻¹와 절대 크기가 약 0.54% 차이로 가깝다.

하지만 ED의 기준 v≈1.07059 nm⁻¹를 b*로 읽으면 XRD b*와 약 4.42% 차이가 있다. ED의 상대 δ≈0.28823과 XRD δ≈0.27826이 다른 이유를 아직 확정하지 않았다. 절대 q의 유사성만으로 같은 3D modulation을 증명할 수 없다. ED 축/setting, calibration, 반사 대응, domain 및 국소 격자 차이를 확인해야 한다. 원래 metadata calibration을 변경하지 않았다.

## 경쟁 가설 및 과적합 점검

1. **동일 hkl + metric 변화만 허용:** B의 다섯 feature를 (001),(020),(110),(-111),(021)로 모두 대응시켜 reciprocal metric 네 계수를 피팅하면 2θ RMS≈0.712°이다. 이 특정 대응의 metric-only 설명은 맞지 않는다. 다른 hkl 대응까지 배제한 것은 아니다.
2. **A의 조건부 metric + 고정 q:** B의 세 추가 feature RMS는 δ=1/4에서 약 0.221°, 1/3에서 0.367°, 1/6에서 0.719°. 1/4가 이 비교군에서는 가장 가깝지만 정밀하게 맞지는 않는다.
3. **A metric + 자유 q + 전체 길이 scale:** B anchor 두 개까지 포함하면 δ≈0.2810, scale≈1.00036에서 RMS≈0.0906°. 일부 feature에는 약 0.17°의 잔차가 남는다.
4. **B metric 네 계수 + δ:** 다섯 위치를 정확히 보간할 수 있다. λ=1.5406 Å라면 대수적 cell은 a≈5.15435, b≈8.94542, c≈5.14372 Å, β≈108.804°다. **5개 관측에 5개 자유도를 쓴 결과로 residual DOF=0이며 refined structure 또는 검증된 격자가 아니다.**

마지막 가설의 별도 예측은 (130)≈35.2462°, (-201)≈35.4167°, (200)≈36.8117°, (002)≈36.8905°다. 이들은 추가로 검증할 fingerprint이며 관측 피크의 hkl 확정값이 아니다.

## 35° 부근의 overlap

A와 B의 35° cluster는 single-envelope 중심이 fitting window에 따라 불안정하고 잔차가 크다. 동일 구간 34.85–35.75°에서 두 envelope를 사용하는 경험적 모델은 잔차 RMS를 A에서 약 116.6→25.1, B에서 218.1→36.2로 줄였다.

- A 두 envelope 중심: 약 35.2013°, 35.3670°.
- B 두 envelope 중심: 약 35.1299°, 35.2803°.
- conditional Gaussian-residual AICc 차이(two−one)는 A≈−125.8, B≈−149.3.

이는 비대칭/겹침을 단일 중심으로 요약하면 위험하다는 증거다. 두 성분이 두 상, 두 hkl, Kα doublet 중 무엇인지는 확정하지 않았다. 처리 이력과 peak-shape가 미확인인 상태에서 AICc로 실제 component 개수를 확정하지 않는다. 따라서 이 feature를 단일 hkl의 정밀 이동값으로 사용하지 않는다. 앞의 free-metric 모델도 이 cluster 전체와의 독립 정합성을 아직 확보하지 못했다.

## 구현·검증·다음 단계

- Workbook reader, local envelope fitting, window/quality flags, wavelength-independent δ, conditional metric hypotheses, overlap diagnostic를 구현했다.
- 자동 발견용 peak detection에만 7-point Savitzky–Golay를 사용했으며 실제 fitting에는 원본 수치를 사용했다.
- 수치 reader, peak center 복원, 파장 불변 δ, overlap 복원 등을 포함해 총 20개 테스트가 통과했다.
- 현재 XRD/ED gate는 not evaluated이며 atomic-structure ranking 및 STEM simulation을 수행하지 않았다.
- 화학적 LFNSO parent CIF와 기존 모델 생성 코드가 있어야 occupancy/displacement/stacking 후보를 원자 수준에서 생성하고 강도·소멸 조건까지 비교할 수 있다. Li₂MnO₃ reference를 임의로 Fe/Ni/Sb 구조로 치환하지 않았다.

현재 실행 결과: `03_structure_model_search/raw_xrd_03/`. `configs/raw_xrd.json`에 입력·파장 가정·피크 window를 저장했다. 재실행은 output을 새 cycle 이름으로 바꾸고 아래 명령을 사용한다.

```bash
MPLCONFIGDIR=/workspace/lfnso-cache/matplotlib PYTHONPATH=src /workspace/lfnso-venv/bin/python -m lfnso.raw_xrd configs/raw_xrd.json
```
