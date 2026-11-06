# MASTER PROMPT — LFNSO21-B Crystallographic Structure Discovery & Multi-Modal Validation

## 0. Project Mission

본 프로젝트의 최종 목표는 **LFNSO21-B의 실험 XRD, ED/SAED, FFT 및 atomic-resolution STEM 결과를 동시에 설명할 수 있는 crystallographic structure model을 체계적으로 탐색하고 검증하는 계산 파이프라인을 구축하는 것**이다.

단일 CIF를 임의로 수정해 실험 영상과 비슷하게 만드는 방식이 아니라, 여러 구조 가설에서 다수의 candidate CIF를 생성하고 실제 실험 데이터와 정량적으로 비교해야 한다.

핵심 전략은 다음과 같다.

**Experimental constraints → Model generation → XRD screening → ED/FFT screening → Model ranking → STEM validation**

특히 계산 비용이 큰 STEM multislice simulation은 XRD와 ED/FFT 검증을 통과한 상위 후보에 대해서만 수행한다.

최종적으로 얻고 싶은 것은 다음과 같다.

1. 실험 결과를 가장 일관성 있게 설명하는 여러 개의 구조 모델과 CIF
2. 후보 구조의 결정학적·화학적 타당성 평가
3. 실험 결과와 구조 모델의 정량적 상관관계
4. 각 구조 모델의 성공·실패 원인 및 구별 가능한 실험적 특징
5. 새로운 자료가 추가되었을 때 탐색과 순위 평가를 반복할 수 있는 reusable codebase

**중요: 이는 단순한 구조 시각화 작업이나 기존 CIF 생성 작업의 반복이 아니다. 구조 후보 탐색·시뮬레이션·실험 대조·모델 선택을 자동화하는 연구용 소프트웨어 개발 작업이다.**

---

# 1. Existing Work and Project Context

## 1.1 기존 작업을 우선 확인

사용자는 이전 Work 프로젝트인 "Develop crystallographic model"과 현재 LFNSO21-B interpretation 프로젝트에서 구조 분석과 CIF 모델링을 반복적으로 수행했다.

Codex는 새로운 구조를 처음부터 임의로 만들기 전에 기존 작업을 확인해야 한다.

현재 로컬 프로젝트 경로:

`C:\Users\user\Desktop\연구\Coding\Codex\LFNSO21-B interpretation`

최근 FFT 분석 결과:

`02_DM_DP_FFT_예비분석\2026-11-05_110_FFT_추가반사_공간분포\report.html`

관련 Work 프로젝트와 로컬 파일에서 기존 분석 결과, CIF, 모델링 코드, 데이터 처리 방법, report 및 simulation output을 탐색한다.

기존에 생성한 것으로 알려진 구조 모델은 다음과 같다.

- `LFNSO21_model0_ideal_ordered.cif`
- `LFNSO21_model1_stack_AABB_registry.cif`
- `LFNSO21_model2_6b_antisite_modulation.cif`
- `LFNSO21_model3_coupled_6b_AABB.cif`
- `LFNSO21_4b_reference_exact_repeat_no_modulation.cif`
- `LFNSO21_4b_q025_TM_antisite_stripe.cif`
- `LFNSO21_4b_q025_displacive_no_stacking_fault.cif`
- `LFNSO21_4b4c_AABB_b6_registry_only.cif`
- `LFNSO21_4b4c_q025_occ75_AABB_b6_coupled.cif`
- 그 외 4b, 6b, stacking fault, occupancy, coupled 모델

기존 CIF는 모두 구조 검증 전의 hypothesis models로 취급한다.

기존 코드의 재사용 가능성을 검토하고, 동일 기능은 중복 구현하지 않는다.

**원본 데이터와 기존 모델은 절대 덮어쓰거나 삭제하지 않는다.**

사용 가능한 파일을 확인한 후 다음을 포함하는 data/model inventory를 작성한다.

- 파일 경로와 형식
- 시료 및 zone axis
- 실험 조건과 calibration
- 데이터 생성 날짜 및 출처
- 기존 구조 모델의 lattice parameters, 조성, site occupancy
- 기존 simulation과 실제 실험의 정합성
- 아직 검증하지 않은 요소

현재 대화 또는 Work 자료에 없는 원자 좌표, XRD peak position, calibration은 추측하여 입력하지 않는다.

---

# 2. Experimental Observations and Structural Constraints

## 2.1 Parent structure

기본 출발점은 Li₂MnO₃-type monoclinic C2/m 기반의 LFNSO21 구조다.

기존 모델링에서 사용했던 lattice parameters는 다음과 같다.

- a = 5.1706 Å
- b = 8.9382 Å
- c = 5.1635 Å
- α = γ = 90°
- β = 109.49°

Nominal composition은 기존 CIF에 사용된 Li₄Fe₂/₃Ni₁/₃SbO₆를 우선 확인한다.

단, 이 수치는 시작점이지 LFNSO21-B에서 확정된 refined structure가 아니다.

참조 CIF의 space group, origin, Wyckoff site, 원소별 occupancy 및 charge neutrality를 확인한다.

## 2.2 XRD constraints

Reference LFSO/LFNSO 또는 LFNSO21-A와 비교했을 때 LFNSO21-B는 다음 특징을 보인다.

- 기존 reference 대비 (110), (111), (021) 등으로 해석된 reflection 위치의 변화
- 기존 A-ordering과 다른 B-ordering superstructure 특징
- B-ordering 관련 reflection이 sharp하고 상대적으로 우세함
- 일부 반사에서는 기존 A-ordering 관련 shoulder 또는 별도 peak 가능성

**반드시 검증할 핵심 질문:**

관찰된 peak shift가 실제 lattice metric 변화인가, 아니면 기존 parent reflection과 다른 superstructure reflection의 출현 및 intensity redistribution인가?

이를 처음부터 결정하지 말고 경쟁 가설로 비교한다.

Bragg peak의 higher 2θ shift는 동일 hkl 반사라는 조건에서 d-spacing contraction을 의미한다.

반면 서로 다른 ordering reflection 또는 중첩된 상의 peak라면 단순 lattice contraction으로 해석해서는 안 된다.

특히 C2/m의 a, b, c, β가 함께 XRD reflection에 영향을 준다는 점을 반영한다.

## 2.3 ED/SAED constraints

**[001] zone axis**

- 기존 parent C2/m pattern과 다른 additional reflection 존재
- (020), (110) 주변에서 satellite-like reflection 관찰
- 기존 자료에서 characteristic separation 약 0.3 nm⁻¹ 수준
- b* 방향 fractional modulation 가능성
- 정확한 q가 1/3b*, 1/4b* 또는 다른 값인지는 미확정

**[101] zone axis**

- (020) 관련 reflection 부근에서 유사한 additional reflection 관찰
- [001]과 동일한 3D modulation vector로 설명 가능한지 검증 필요

**[110] zone axis**

- 기존 LFSO/LFNSO보다 두 줄의 diffuse streak가 감소하는 경향
- 비교적 sharp한 spot 형태
- 높은 reciprocal-space frequency에서도 기존 parent model과 다른 additional FFT components
- Line-like HAADF contrast와의 공간적 연관성 검증 필요

**기존 [100]/[010] 및 기타 회절자료**

- 해당 자료가 존재하면 독립적인 structural constraints로 포함
- [010]에서 보고된 reflection doublet과 기존 C2/m metric의 정합성 검토
- 단, 같은 raw data에서 생성된 중복 분석 결과를 독립적인 실험으로 취급하지 말 것

## 2.4 Atomic-resolution STEM constraints

**[001] projection**

- 기존 honeycomb-type ordering과 다른 projected atomic-column arrangement
- 약 3 strong–3 weak columns 또는 4–5 atomic-column motif처럼 보이는 특징
- 일정한 spatial modulation 가능성
- 하지만 원소 점유, 원자 위치 변위, 여러 slab의 projection overlap을 구별해야 함

**[110] projection**

- 기존 dumbbell ordering과 다른 line-like contrast
- 동일 입자의 서로 다른 영역에서 dot-like / line-like contrast 가능성
- TM slab의 layered framework는 비교적 연속적으로 유지
- 일부 column이 slab 내에서 미세하게 위아래로 displacement된 것처럼 관찰됨
- 해당 displacement가 실제 구조적인지, STEM scan distortion인지 검증 필요

현재 관찰된 특징을 강제로 한 가지 구조의 증거라고 가정하지 않는다.

---

# 3. Research Hypotheses to Test

다음 가설을 모두 모델링하되, 가능한 한 구조적 자유도를 독립적으로 변화시켜 각 효과를 분리한다.

### H0 — Average metric distortion only

기존 C2/m honeycomb arrangement와 stacking을 유지하고 lattice parameters만 변화시킨다.

목적: XRD peak shift가 평균 격자 변화만으로 설명되는지 확인.

### H1 — Reduced random stacking disorder

기존 honeycomb ordering을 유지하되 stacking fault density와 stacking correlation length를 변화시킨다.

목적: Diffuse streak 감소 및 reflection sharpening 설명.

### H2 — Correlated/ordered stacking registry

개별 TM slab의 in-plane ordering을 유지하고 slab 간 lateral registry 및 stacking sequence를 변화시킨다.

목적: [110] line-like contrast, stacking-sensitive reflections 및 XRD B-ordering 변화 설명.

### H3 — In-plane occupational modulation

TM slab 내 Li/Fe/Ni/Sb occupancy 또는 cation ordering을 주기적으로 변화시킨다.

목적: [001]/[101] fractional-b* reflection 및 projected atomic-column contrast 설명.

### H4 — Periodic atomic displacement

기존 occupancy를 유지하면서 원자 위치의 주기적인 변위를 도입한다.

목적: ED satellite와 [110] column displacement 설명.

### H5 — Coupled occupational and displacive modulation

Periodic occupancy modulation에 따른 TM/O sublattice relaxation을 추가한다.

목적: 추가 반사의 intensity dependence와 실공간 왜곡을 동시에 설명.

### H6 — Coupled in-plane ordering and correlated stacking

In-plane ordering과 layer-dependent ordering phase 및 stacking registry를 함께 변화시킨다.

현재 가장 중요하게 검증하고 싶은 구조 가설이다.

### H7 — Coexisting structural domains

동일 입자 또는 powder sample 내부에 서로 다른 ordering domains, twin variants 또는 coherent/incoherent structural variants가 공존할 가능성.

### H8 — Alternative polytypism / intergrowth

필요한 경우 O1/O3-like stacking 또는 다른 layered polytype을 검토한다.

단, H0–H7의 상대적으로 단순한 모델이 실패했을 때 우선순위를 높인다.

**중요: 가장 복잡한 H6/H8을 처음부터 정답으로 선택하지 않는다. 단순 모델을 대조군으로 유지하면서 추가 자유도가 실제로 필요한지 증명한다.**

---

# 4. Automated Crystallographic Model Generation

## 4.1 General modeling rules

기존 CIF를 기반으로 symmetry를 적절히 해제하고 supercell을 생성한다.

실제 site occupancy를 변경하는 과정에서는 필요할 때 P1 representation을 이용하되, 계산 편의를 위한 P1 설정이 실제 결정의 space group 확정으로 해석되지 않도록 한다.

각 모델에 다음 정보를 기록한다.

- Parent CIF
- Lattice transformation matrix
- Site modification history
- Total chemical composition
- Occupancy and charge constraints
- Atomic displacement vectors
- Stacking sequence
- Physical superperiodicity
- Candidate space-group symmetry
- Model generation seed and parameters

## 4.2 In-plane supercell candidates

최소한 다음을 포함한다.

- 1b reference
- 3b supercell
- 4b supercell
- 6b supercell
- 3b/4b와 비교 가능한 12b approximant
- 필요시 10b 또는 다른 rational approximant

기존 b = 8.9382 Å 기준,

- 3b: q = 0.373 nm⁻¹
- 4b: q = 0.280 nm⁻¹
- 6b: q = 0.186 nm⁻¹

이들은 검증할 후보이며 실험적으로 확정된 값이 아니다.

관찰된 spot separation이 fundamental q인지, 2q harmonic인지, domain splitting인지부터 판별한다.

Large supercell을 처음부터 모두 계산하지 말고, 작은 모델에서 설명되지 않는 peak position을 다룰 때 단계적으로 확장한다.

## 4.3 Occupational modulation parameters

서로 다른 site 유형을 구분하여 변경한다.

**TM slab 내 site**

- Li-rich site
- Fe/Ni mixed site
- Sb site
- 필요시 Fe/Ni ordering variation

**Interlayer Li site**

- Interlayer occupancy control
- TM antisite control
- Layer-dependent occupation control

실제 Wyckoff position은 기존 CIF를 읽어 확인해야 하며 이름만으로 site를 추정하지 않는다.

다음 변수를 탐색한다.

- Li/Fe exchange fraction
- Li/Ni exchange fraction
- Fe/Ni ordering ratio
- Li/(Fe,Ni) antisite concentration
- Sb site occupancy 변동을 포함하는 제한된 대조 모델
- Periodic occupational stripe
- Ordering-domain phase
- Layer-dependent composition
- Symmetry-related ordering variants

예시 occupational exchange fractions:

`0, 0.10, 0.25, 0.50, 0.75, 1.00`

이 값들은 탐색용 초기 범위이며, 실제로 불가능한 화학 상태나 조성은 구조 검증 단계에서 탈락시킨다.

단순히 모든 Li site를 TM으로 치환하는 모델은 별도 control로만 사용한다.

**Stoichiometry-conserving site exchange를 기본 모델로 우선 사용한다.**

Partially occupied site 모델과 실제 정수 원자 배열을 갖는 ordered approximant를 구분한다. 부분 점유의 평균 구조와 하나의 구체적인 원자 배열을 동일하게 취급하지 않는다.

## 4.4 Atomic displacement modulation

다음 displacement mode를 별도로 구현한다.

- Longitudinal displacement along modulation direction
- Transverse displacement
- Slab-normal displacement
- TM sublattice buckling
- Oxygen sublattice relaxation
- Correlated TM–O displacement
- Layer-dependent displacement phase

예시:

`u(y) = A cos(2πδy/b + φ)`

변조 방향과 원자 displacement direction을 별개의 매개변수로 취급한다.

Displacement amplitude는 낮은 값부터 시작해 점진적으로 탐색한다.

원자 간 비정상적으로 짧은 결합거리, 비현실적인 octahedral geometry, layer overlap 등은 자동 검출한다.

## 4.5 Interlayer stacking modulation

다음 구조를 탐색한다.

- Conventional parent stacking
- AABB-like stacking
- ABAB-like stacking
- AAAB 또는 기타 짧은 주기 배열
- Layer registry shift
- Layer-dependent ordering phase
- Periodic stacking faults
- Random stacking faults
- Correlated stacking disorder
- Rotational twins
- Alternative layered polytypes

기존에 제안된 ±b/6 및 a 방향 translation을 포함하되 이것을 유일한 stacking operation으로 가정하지 않는다.

Fractional translation은 실제 lattice-vector convention에 따라 정확히 변환한다.

각 slab의 translation 시 O sublattice 및 edge-sharing octahedral connectivity까지 고려한다. TM cation만 임의로 이동시켜 화학적으로 불가능한 구조를 생성하지 않도록 한다.

Disorder 모델은 단순한 완전 주기적 supercell만으로 대체하지 않는다. 필요하면 finite stacking sequence 또는 ensemble averaging을 도입한다.

## 4.6 Combined parameter space

다음 세 종류의 자유도를 독립적으로 구성한다.

- Occupational parameters: O
- Displacement parameters: D
- Stacking parameters: S

모델 집합:

`M(O,D,S)`

최소 비교 그룹:

- M(0,0,0)
- M(O,0,0)
- M(0,D,0)
- M(0,0,S)
- M(O,D,0)
- M(O,0,S)
- M(0,D,S)
- M(O,D,S)

단, 모든 조합을 무작위로 전수조사하지 않는다.

**Hierarchical search + coarse-to-fine optimization을 구현한다.**

---

# 5. Stage 1 — XRD-Based Screening

## 5.1 Average parent lattice refinement

사용 가능한 experimental XRD raw data를 확인한다.

기본 C2/m indexing을 기준으로 fundamental reflections를 선정하고,

`a, b, c, β`

를 refinement한다.

필요에 따라 zero shift, instrumental broadening, Kα splitting, preferred orientation 등의 영향을 고려한다.

관찰된 (110), (111), (021)의 shift를 무조건 lattice contraction이라고 가정하지 않는다.

특히 A/B superstructure reflection과 fundamental reflection을 구별해야 한다.

**필수 비교:**

- Uniform lattice distortion model
- New superstructure reflection model
- Mixed-domain / multiphase model
- Stacking-related broadening or intensity redistribution model

반사 하나의 이동을 다른 hkl로 잘못 대응시키는 peak misassignment를 자동 검토한다.

Fundamental reflections가 충분하지 않으면 4개의 monoclinic metric parameter를 모두 유일하게 refinement할 수 있다고 주장하지 않는다.

## 5.2 Powder diffraction simulation

가능한 경우 pymatgen XRDCalculator 등을 이용하여 후보 구조의 powder diffraction pattern을 계산한다.

다음 항목을 평가한다.

- Fundamental peak positions
- Superstructure peak positions
- A/B ordering peak separation
- Integrated intensity ratios
- Peak splitting
- Missing or forbidden reflections
- Unexpected strong reflections

Ideal periodic CIF 기반 XRD가 stacking disorder에 의한 diffuse/broadening을 설명하지 못하는 경우, finite stacking ensemble 또는 disorder-aware diffraction model을 별도로 사용한다.

## 5.3 XRD scoring

실험 자료의 품질과 오차를 반영하여 평가한다.

예시 평가항목:

- Fundamental peak position error
- Superstructure peak position error
- Relative intensity mismatch
- Profile residual
- Unobserved strong peak penalty
- Model complexity penalty

측정되지 않은 peak나 detection limit 아래의 peak를 임의로 0 intensity라고 가정하지 않는다.

실험 XRD가 없는 경우 가상의 데이터를 생성하지 않고 해당 constraint를 unavailable로 표시한다.

**Gate 1:** XRD와 구조적·화학적으로 양립할 수 없는 모델을 탈락시킨다.

단, XRD에서 불확실하거나 domain averaging으로 영향을 받는 항목은 확정적인 exclusion criterion으로 사용하지 않는다.

---

# 6. Stage 2 — ED / SAED / FFT Screening

## 6.1 Accurate crystallographic indexing

Experimental ED를 기존 C2/m parent reciprocal basis로 indexing한다.

Supercell simulation에서는 index 변환을 정확히 적용한다.

예를 들어 3b supercell에서는:

- Parent (020) → Supercell (060)
- Parent (110) → Supercell (130)
- Parent [110] ZA → Supercell [310] ZA

4b에서는:

- Parent (020) → Supercell (080)
- Parent (110) → Supercell (140)
- Parent [110] ZA → Supercell [410] ZA

실제 reciprocal vector가 같으면 index가 달라도 동일한 reflection으로 매칭한다.

별도 notation으로 parent G와 modulation q를 유지한다.

## 6.2 Reciprocal-space peak extraction

원본 DM3/DM4, TIFF, FFT data에서 peak coordinate를 정량적으로 추출한다.

필수 기능:

- Reciprocal calibration
- Direct-beam center fitting
- Parent lattice fitting
- Satellite peak detection
- Peak broadening / FWHM extraction
- Local background subtraction
- Peak uncertainty estimation
- Reciprocal vector transformation
- Symmetry-related peak comparison

Raw peak와 FFT에서 생성된 peak를 구별한다.

특히 실제 TEM diffraction intensity와 HAADF FFT intensity는 동일한 물리량이 아니므로 같은 scoring rule을 직접 적용하지 않는다.

## 6.3 Global modulation-vector analysis

[001], [101], [110], 가능하다면 [100]/[010]에서 추출한 추가 반사를 하나의 3D reciprocal model에 매칭한다.

후보:

- q ≈ 1/3 b*
- q ≈ 1/4 b*
- q ≈ 1/6 b*
- Mixed q-vector
- q와 2q harmonic
- Incommensurate modulation
- Twinned/domain-overlap patterns
- Stacking-sensitive reflections

단순히 q의 절대값이 같다는 이유로 동일한 modulation이라고 판단하지 않는다.

정확한 vector direction, reciprocal basis, zone-law condition 및 parent indexing을 확인한다.

[001]/[101]에서 관찰한 순수 b* modulation이 이상적인 [110] ZOLZ에서 직접 관찰되는 반사인지 검증하고, 그렇지 않은 경우 [110] extra reflection을 별도 구조적 constraint로 다룬다.

## 6.4 Diffraction simulation

모든 모델에 대해 먼저 빠른 kinematic structure-factor 계산을 수행한다.

대량 후보 검색 단계에서는 다음을 우선 비교한다.

- Reflection positions
- Expected extinctions
- Satellite positions
- Relative structure-factor intensities
- Harmonic generation
- Multiple equivalent reciprocal lattice points

3b/4b supercell에서 단순 반복만 수행한 대조 모델에는 인위적인 superlattice intensity가 나타나지 않아야 한다.

해당 테스트는 필수 unit test로 구현한다.

실험과 기본적인 정합성이 있는 모델에 대해서만 더 정교한 electron diffraction simulation을 수행한다.

Electron diffraction intensity가 dynamical scattering, thickness, orientation, excitation error의 영향을 받는다는 것을 반영한다.

JEMS와 비교할 수 있는 reflection tables 및 simulated patterns도 출력한다.

## 6.5 FFT-specific analysis

[110] 고차 FFT에서 관찰한 추가 spot을 다음과 같이 분류한다.

- Parent reflection
- Harmonic of parent structure
- New superstructure reflection
- Satellite reflection
- Stacking-related reciprocal intensity
- Twin/domain contribution
- Scan artifact 또는 image processing artifact

FFT/IFFT를 사용하여 다음을 비교한다.

- Parent-only reconstruction
- Additional-reflection-only reconstruction
- Parent + Additional reconstruction
- Residual signal

FFT mask는 Friedel symmetry와 conjugate pairing을 고려한다.

Mask 크기와 위치에 따라 생기는 인공적인 stripe 또는 dot-like contrast를 정량적으로 평가한다.

Local FFT 결과를 해석할 때에는 field-of-view에 의한 reciprocal resolution과 window leakage를 고려한다.

**Gate 2:** XRD 및 여러 zone의 reciprocal-space 결과를 동시에 설명하는 모델만 최종 STEM 후보로 선정한다.

---

# 7. Structural Search and Optimization Engine

## 7.1 Search strategy

수많은 supercell 후보를 전부 expensive simulation으로 평가하지 않는다.

다음 방식의 탐색 시스템을 구현한다.

**Phase A — Screening**

- Symmetry-distinct configurations 생성
- Chemical/geometric validity 검사
- Fast structure factor
- Approximate XRD/ED 비교

**Phase B — Optimization**

- 살아남은 candidate에 대해 occupancy fraction, modulation amplitude, phase, stacking sequence refinement
- 필요하면 Bayesian optimization, evolutionary optimization 또는 multi-objective optimization 적용
- 동일 구조에 대한 중복 계산 방지

**Phase C — High-fidelity validation**

- 상위 모델에 대한 상세 XRD comparison
- Multi-zone electron diffraction simulation
- 필요시 disorder ensemble
- 이후 STEM simulation

초기 탐색의 모델 수와 계산 예산은 사용 가능한 하드웨어를 조사한 후 설정한다.

GPU가 없어도 low-cost screening 단계는 작동해야 한다.

## 7.2 Model ranking

하나의 가중 합산 score만으로 최적 구조를 결정하지 않는다.

다음 요소를 별도로 평가한다.

- XRD consistency
- ED consistency
- FFT consistency
- Physical validity
- Model complexity
- Parameter uncertainty

Pareto ranking 또는 유사한 multi-objective comparison을 구현한다.

실험 데이터의 품질이 불균등할 수 있으므로 각 항목의 uncertainty를 고려하고, 근거 없는 가중치로 구조 순위를 임의 조정하지 않는다.

동일한 결과를 설명하는 여러 구조는 crystallographically indistinguishable한 후보군으로 남긴다.

## 7.3 Overfitting prevention

가능하다면 사용한 반사군과 검증에 사용할 반사군을 분리한다.

예를 들어 일부 XRD reflection과 [001] ED를 fitting한 뒤, 다른 hkl 및 [101] ED를 hold-out validation으로 사용할 수 있다.

특정 후보가 하나의 zone에만 잘 맞고 다른 zone에 실패하면 그 한계를 반드시 기록한다.

다음 negative controls를 필수적으로 포함한다.

- Exact-repeat supercell
- Ideal ordered parent
- Stacking-only model
- Occupancy-only model
- Displacement-only model
- Randomized occupancy control
- Domain-overlap control

---

# 8. Stage 3 — STEM Validation (Shortlisted Models Only)

**STEM simulation은 Gate 1과 Gate 2를 통과한 후보만 진행한다.**

검증 방향:

- [001]
- [110]
- 필요시 [100]

가능하다면 abTEM multislice/PRISM 또는 기존에 사용하던 적절한 STEM simulation package를 이용한다.

## 8.1 [001] validation

다음의 projected characteristics를 비교한다.

- Honeycomb-type ordering
- 3 strong–3 weak column motif
- 4–5 atomic-column periodicity
- Column intensity profile
- Local superstructure modulation
- Projection of multiple TM slabs

## 8.2 [110] validation

다음 결과를 평가한다.

- Dumbbell-like vs line-like contrast
- Intensity of nominal weak columns
- TM slab continuity
- Column lateral displacement
- Slab buckling
- Thickness-dependent projection
- Different ordering domains

가능한 경우 시료 두께, orientation deviation, probe parameters, detector angles, defocus/source size 등의 변수를 분리하여 평가한다.

실험 조건이 없는 경우 적절한 범위의 sensitivity analysis로만 사용하며 특정 조건을 실제 측정값처럼 기재하지 않는다.

## 8.3 Quantitative comparison

가능한 지표:

- Atomic column position error
- Intensity profile similarity
- Projected periodicity
- Local displacement statistics
- FFT peak position and relative intensity
- Structural motif similarity

Pixel-wise image correlation만으로 구조의 우수성을 판단하지 않는다.

Scan distortion, sample tilt, thickness 및 projection overlap을 별도 nuisance parameters로 취급한다.

STEM 단계에서 구조 후보가 실패하면 해당 모델을 기록하고 재평가한다. 고비용 STEM 계산을 모든 초기 후보에 반복하지 않는다.

---

# 9. Creative Model Extensions

기본 탐색만으로 결과가 설명되지 않는 경우 다음 아이디어를 추가한다.

## 9.1 Ordering-phase-coupled stacking

모든 TM slab이 동일한 occupational modulation을 갖는 대신, layer에 따라 modulation phase가 달라지는 모델.

예시:

`p(y,n) = p0 + A cos(2πδy/b + φn)`

여기서 φn은 층마다 변화할 수 있다.

AABB stacking뿐 아니라 phase-locked 또는 antiphase ordering을 포함한다.

이 모델이 diffraction intensity와 [001] projected contrast에 어떤 영향을 주는지 계산한다.

## 9.2 Interlayer projection vs genuine in-plane change

[001]에서 동일한 3 strong–3 weak contrast를 만드는 서로 다른 구조를 생성한다.

- TM slab 내 occupational rearrangement
- Interlayer TM antisite
- Stacking registry variation
- Multiple ordering domains

[110] 구조와 ED를 이용해 구별한다.

## 9.3 Correlated buckling

TM slab 내부의 periodic cation ordering이 oxygen octahedral distortion 또는 TM column buckling과 연동되는 모델.

Occupancy-only, displacement-only, coupled models의 diffraction fingerprints를 비교한다.

## 9.4 Domain-population model

XRD는 bulk-averaged data인 반면 STEM은 국소 구조를 반영한다.

따라서 모든 입자에 동일한 구조가 존재해야 한다고 가정하지 않는다.

필요하면 서로 다른 structural-domain fraction을 도입해 powder pattern을 계산하고, ED의 국소 구조와 일치하는지 확인한다.

단, 여러 모델의 혼합 비율을 자유롭게 늘리는 방식으로 과적합하지 않도록 complexity penalty를 적용한다.

## 9.5 Diffraction fingerprint library

각 후보 구조마다 다음 특징을 자동 기록한다.

- Strongest distinguishing reflection
- Missing reflection
- Predicted q vector
- [001]/[101]/[110] characteristic peaks
- A/B-ordering XRD intensity ratio
- Sensitivity to occupancy/displacement/stacking

최종적으로 실험 결과로 모델을 구분하는 데 유용한 reflection과 zone axis를 자동 제안한다.

---

# 10. Code Architecture and Reproducibility

Python 기반으로 modular package를 구성한다.

권장 기능 단위:

- `data_inventory`
- `data_io`
- `crystal_geometry`
- `parent_refinement`
- `model_generator`
- `chemical_validation`
- `xrd_simulation`
- `ed_simulation`
- `fft_analysis`
- `search_optimizer`
- `model_ranking`
- `stem_validation`
- `report_generator`

사용할 수 있는 라이브러리를 먼저 확인한다.

예시:

- pymatgen
- ASE
- spglib
- gemmi
- NumPy / SciPy
- pandas
- NCEMpy 또는 HyperSpy
- matplotlib / Plotly
- abTEM

버전과 API compatibility를 확인한 뒤 채택한다.

필요한 기능을 지원하지 않는 라이브러리를 억지로 적용하지 않는다.

상용 JEMS의 기능은 자동 실행 가능 여부를 확인하고, 접근할 수 없으면 비교용 reflection table 및 CIF export를 생성한다.

## 10.1 Required reproducibility features

- Configuration-driven execution
- Deterministic random seed
- Input data manifest
- Model parameter logging
- Structural validity tests
- Calculation caching
- Checkpoint and resume
- Parallelization where safe
- Intermediate results export
- No modification of raw data

기존 파일을 사용해 계산할 경우 원본 경로와 모델 provenance를 항상 기록한다.

## 10.2 Output organization

새로운 프로젝트 하위 디렉터리를 생성한다.

예시:

`03_structure_model_search/`

내부에는 다음을 구분한다.

- `inputs_manifest/`
- `reference_structures/`
- `candidate_cifs/`
- `xrd_results/`
- `ed_fft_results/`
- `rankings/`
- `stem_shortlist/`
- `reports/`
- `configs/`
- `logs/`

생성 파일과 기존 파일이 섞이지 않도록 한다.

---

# 11. Required Final Reports

각 탐색 cycle에서 자동으로 HTML report와 machine-readable 결과를 생성한다.

최소 내용:

### A. Data quality and constraints

- 어떤 experimental data를 사용했는가?
- 어떤 데이터가 부족한가?
- 각 데이터의 uncertainty는 얼마인가?

### B. Model search summary

- 생성한 모델 수
- 화학적·결정학적 이유로 탈락한 모델 수
- XRD/ED gate 통과 모델 수
- 계산 시간 및 resource usage

### C. Ranked candidates

각 상위 후보에 대해:

- CIF
- Structure parameters
- Chemical composition
- Space group / supercell
- Occupational ordering
- Stacking sequence
- Displacement modulation
- XRD/ED agreement
- 실패한 experimental constraints

### D. Direct experimental comparison

- Experimental vs simulated XRD
- Experimental vs simulated [001] ED
- Experimental vs simulated [101] ED
- [110] additional reflection comparison
- Satellite positions 및 intensity ratios
- Predicted vs experimental FFT features

### E. Scientific interpretation

다음을 명확히 보고한다.

1. 어떤 구조적 변화가 실험 데이터를 설명하는 데 필수적인가?
2. In-plane ordering 변화 없이 stacking만으로 설명할 수 있는가?
3. Stacking 변화 없이 in-plane modulation만으로 충분한가?
4. Atomic displacement를 도입해야 하는가?
5. 서로 구별되지 않는 구조 후보는 무엇인가?
6. 최적 후보를 반증할 수 있는 실험 결과는 무엇인가?

### F. Recommendations

모델 간 차이가 가장 크게 나타나는 reflection 또는 zone axis를 찾아, 추가 실험이 필요할 경우 효율적인 검증 방향을 제안한다.

단순히 "추가 TEM 측정이 필요하다"는 결론을 반복하지 말고 구체적인 discriminating signatures를 제시한다.

---

# 12. Implementation Instructions — Start Actual Work

**코딩 설계만 작성하고 종료하지 말 것. 실제 실행 가능한 프로그램을 구현하고 기존 실험 데이터로 첫 번째 탐색을 수행한다.**

다음 순서로 작업한다.

### Milestone 1 — Existing project audit

기존 파일, CIF, 분석 코드 및 report를 확인한다.

특히 기존 4b/6b occupational, displacive, AABB stacking 모델의 구현 방법과 검증 결과를 파악한다.

누락된 입력 데이터는 명확히 표시한다.

### Milestone 2 — Unified crystallographic indexing

기존 parent C2/m과 다양한 supercell의 reciprocal transformation을 구현한다.

기존 XRD/ED reflection을 공통 reciprocal-space convention으로 비교한다.

자동 테스트를 포함한다.

### Milestone 3 — Automated model generation

Reference, occupancy-only, displacement-only, stacking-only 및 coupled structure를 자동 생성한다.

CIF export와 조성·기하학적 검증을 수행한다.

### Milestone 4 — First screening

사용 가능한 XRD와 ED/FFT 데이터를 바탕으로 초기 구조 검색을 실행한다.

기본 대조 모델에서 시작해 occupancy와 stacking parameters를 단계적으로 확장한다.

### Milestone 5 — Ranking and scientific report

상위 모델과 실패한 모델의 주요 차이를 정량적으로 보고한다.

사용자의 확인 없이 특정 구조가 확정되었다고 선언하지 않는다.

### Milestone 6 — Advanced validation

XRD 및 ED 검증을 통과한 후보만 대상으로 추가적인 diffraction refinement와 STEM multislice validation을 실행한다.

---

# 13. Critical Research Principles

반드시 다음 원칙을 지킨다.

1. **Do not assume a unique crystal structure.** 서로 다른 구조가 유사한 projected contrast와 diffraction을 만들 수 있음을 고려한다.

2. **Do not equate a supercell with physical modulation.** 동일 구조를 단순 반복한 supercell에서는 새로운 superlattice intensity가 생기지 않아야 한다.

3. **Do not confuse reciprocal indexing changes with physical peak shifts.** Parent와 supercell의 reflection index는 서로 다르더라도 실제 reciprocal vector가 같을 수 있다.

4. **Do not infer chemical occupancy directly from HAADF intensity alone.** Thickness, projected overlap, detector response 및 tilt를 고려해야 한다.

5. **Do not attribute every extra FFT spot to a new superstructure.** Harmonics, stacking-related components, twins, projection effects 및 artifacts를 구별한다.

6. **Do not fit only a single experimental observation.** XRD, multi-zone ED 및 STEM의 정합성을 우선한다.

7. **Avoid unphysical modifications.** Chemical stoichiometry, site occupancy, interatomic distances, oxygen coordination 및 slab connectivity를 검증한다.

8. **Prefer the simplest model that explains all reliable constraints.** 추가적인 자유도가 실제 설명력을 개선하는지 검증한다.

9. **Quantify uncertainty and model ambiguity.** 실험적으로 구별되지 않는 모델은 함께 남긴다.

10. **Preserve previous work and research continuity.** 기존 CIF와 Work 결과를 재사용하고, 원본 입력을 훼손하지 않는다.

---

# 14. Final Research Goal

최종적으로 LFNSO21-B에서 다음 구조적 질문에 답하는 것이 목표다.

**Is the observed B-ordering superstructure primarily governed by modified in-plane cation ordering, correlated interlayer stacking, periodic lattice distortion, or a coupled combination of these structural degrees of freedom?**

가능한 구조를 제한 없이 늘리는 것이 목적이 아니다.

실험 데이터를 가장 잘 설명하면서도 화학적·결정학적으로 타당한 소수의 모델을 도출하고, 각 모델의 성공 및 실패 원인을 정량적으로 제시해야 한다.

현재 가장 중요하게 검증할 가설은 다음과 같다.

**A correlated 3D ordering state in which modified in-plane cation arrangements coexist with more coherent interlayer stacking, possibly accompanied by periodic lattice relaxation.**

그러나 이것은 출발 가설일 뿐이다. 계산 결과가 이 가설과 다르면 다른 구조를 선택해야 한다.

최종적으로 사용자가 얻어야 하는 것은 단순한 CIF 파일 여러 개가 아니라, **새로운 experimental XRD/ED/STEM data를 입력할 때마다 구조 후보를 생성·검증·갱신할 수 있는 crystallographic structure discovery framework**이다.

지금부터 기존 프로젝트를 조사하고, 가장 작은 end-to-end pipeline을 먼저 실제로 구현·검증한 다음 단계적으로 모델 탐색 범위를 확장하라. 실제 입력 데이터가 허용하는 범위 내에서 첫 번째 ranked candidate 결과까지 도출하라.