# 🤖 Agent Handoff — quantum-protein-folding

> **Dla:** Claude Opus 5  
> **Branch:** `feat/research-project-2s`  
> **Stan testów:** `312 passed` ✅ (`uv run pytest tests/ -q`)  
> **Data:** 2026-09-09  
> **Repozytorium:** `c:\repos\quantum-protein-folding`

---

## 1. Czym jest ten projekt

**Quantum Protein Folding** to implementacja algorytmu kwantowego składania białek opartego na siatce tetraedrycznej (FCC/diamond lattice). Białko kodowane jest jako ciąg turn-qubitów opisujących kolejne skręty na kratce przestrzennej. Zadaniem jest minimalizacja Hamiltonianu energetycznego przy użyciu algorytmu VQE.

### Model fizyczny

```
H_total = H_backbone + H_backtrack + H_interaction + H_field (opt.) + H_ligand (opt.)
```

- **H_backbone** — penalizuje niedozwolone ruchy wzdłuż łańcucha  
- **H_backtrack** — eliminuje samoprzecięcia łańcucha  
- **H_interaction** — model oddziaływań HP lub Miyazawa-Jernigan (MJ) między aminokwasami  
- **H_field** — zewnętrzne pole energetyczne (Wariant A: globalne przesunięcie skalarne per bead)  
- **H_ligand** — oddziaływanie ligand–białko (Wariant A: scalar shift sumowany po wszystkich pozycjach kratki)

### Kodowanie qubitowe

Każde obracanie łańcucha zajmuje `QUBITS_PER_TURN = 2` qubity (stała tetrahedral FCC). Łańcuch N aminokwasów = `(N-1) * 2` turn-qubitów + ewentualne qubity liganda.

---

## 2. Struktura projektu

```
c:\repos\quantum-protein-folding\
├── src/
│   ├── main.py                       # Entry point: APRLRFY, pełny workflow VQE
│   ├── enums.py                      # InteractionType (HP, MJ), Penalties
│   ├── constants.py                  # QUBITS_PER_TURN=2, energie kontaktów, ścieżki
│   ├── exceptions.py                 # Własne wyjątki (InvalidInteractionTypeError, etc.)
│   │
│   ├── particle/
│   │   ├── external_field.py         # [A1] ExternalField (UNIFORM/NON_UNIFORM)
│   │   └── ligand_bead.py            # [B1] LigandBead (BINARY/UNARY pozycja kwantowa)
│   │
│   ├── interaction/
│   │   ├── interaction.py            # ABC: get_energy(sym_i, sym_j) → float
│   │   ├── hp_interaction.py         # HPInteraction: H-H=-1.0, reszta=0.0
│   │   ├── mj_interaction.py         # MJInteraction: macierz 20x20 aminokwasów
│   │   └── ligand_interaction.py     # [B2] LigandInteraction (HP_LIKE / CUSTOM)
│   │
│   ├── builder/
│   │   └── hamiltonian_builder.py    # [A2,B3] HamiltonianBuilder: sum_hamiltonians()
│   │
│   ├── analysis/
│   │   ├── field_influence_analysis.py  # [A3] FieldInfluenceAnalysis: sweep λ
│   │   └── ligand_analysis.py           # [B4] LigandAnalysis: VQE + binding energy
│   │
│   ├── protein/                      # Protein, Bead (turn-qubit chain)
│   ├── contact/                      # ContactMap — które beady sąsiadują?
│   ├── distance/                     # DistanceMap — odległości na kratce
│   ├── optimization/                 # VQE setup utilities
│   ├── backend/                      # get_sampler() → AerSimulator
│   ├── result/                       # Interpretacja i wizualizacja wyników VQE
│   ├── visualization/                # 2D/3D, GIF konformacji białka
│   ├── utils/
│   │   ├── qubit_utils.py            # build_turn_qubit, build_identity_op, pad_to_n_qubits
│   │   └── setup_utils.py            # Konveniencyjne funkcje dla main.py
│   ├── validation/                   # Walidacja sekwencji aminokwasów
│   ├── logger/                       # get_logger() — standardowy logger
│   └── resources/                    # Pliki danych: hp_interaction_matrix.txt, mj_matrix.txt
│
├── tests/
│   ├── test_external_field.py        # [A1] 41 testów
│   ├── test_hamiltonian_external_field.py  # [A2] 13 testów
│   ├── test_field_influence_analysis.py    # [A3] 35 testów
│   ├── test_ligand_bead.py           # [B1] 49 testów
│   ├── test_ligand_interaction.py    # [B2] 82 testy
│   ├── test_hamiltonian_ligand.py    # [B3] 23 testy
│   ├── test_ligand_analysis.py       # [B4] 55 testów
│   └── test_utils.py                 # 14 testów
│
├── docs/
│   └── field_extension_presentation.ipynb  # Notebook: pełna prezentacja Wariantu A+B
│
├── notebooks/                        # Inne notebooki (eksploracyjne)
├── output/                           # Domyślny katalog wyjściowy (wykresy, GIF-y)
├── pyproject.toml                    # Zależności: qiskit, qiskit-algorithms, scipy, etc.
├── walkthrough.md                    # Szczegółowy opis etapów B1–B4
└── AGENT_HANDOFF.md                  # Ten plik
```

---

## 3. Co zostało zrobione — historia etapów

### Etap A1 — ExternalField ✅

**Plik:** [`src/particle/external_field.py`](file:///c:\repos\quantum-protein-folding\src\particle\external_field.py)

Klasa `ExternalField` reprezentuje zewnętrzne pole energetyczne działające na łańcuch aminokwasów.

```python
from particle.external_field import ExternalField, FieldMode

# Tryb jednorodny: każdy węzeł kratki = ta sama energia
field = ExternalField.uniform(strength=-1.0)
field.get_energy((3,))   # → -1.0

# Tryb niejednorodny: słownik lattice_coords → energy
field2 = ExternalField.non_uniform({(0,): -2.0, (1,): 0.5, (3,): -1.5})
field2.get_energy((1,))  # → 0.5
field2.get_energy((9,))  # → 0.0  (domyślna, nieznany węzeł)

# Nadpisanie węzła
field.set_energy((5,), -3.0)
```

**41 testów** w `tests/test_external_field.py`.

---

### Etap A2 — HamiltonianBuilder + H_field ✅

**Plik:** [`src/builder/hamiltonian_builder.py`](file:///c:\repos\quantum-protein-folding\src\builder\hamiltonian_builder.py)

Rozszerzony `HamiltonianBuilder.__init__` o parametr `external_field: ExternalField | None = None`.

```python
from builder.hamiltonian_builder import HamiltonianBuilder
from particle.external_field import ExternalField

builder = HamiltonianBuilder(
    protein=protein,
    interaction=interaction,
    distance_map=distance_map,
    contact_map=contact_map,
    external_field=ExternalField.uniform(-1.0),   # ← nowy parametr
)
H = builder.sum_hamiltonians()
```

**Ważne — ograniczenie Wariantu A:**  
`H_field = Σ_i E_field((i,)) · I` — skalarne przesunięcie tożsamościowe. Przesuwa energie WSZYSTKICH stanów o tę samą stałą. Optymalna konformacja białka NIE zmienia się przy różnych λ — relatywne różnice energii między stanami są zachowane.

**13 testów** w `tests/test_hamiltonian_external_field.py`.

---

### Etap A3 — FieldInfluenceAnalysis ✅

**Plik:** [`src/analysis/field_influence_analysis.py`](file:///c:\repos\quantum-protein-folding\src\analysis\field_influence_analysis.py)

Klasa orkiestrująca VQE dla zestawu scenariuszy pola zewnętrznego.

```python
from analysis.field_influence_analysis import FieldInfluenceAnalysis
from enums import InteractionType

analysis = FieldInfluenceAnalysis(
    main_chain="APRLRFY",
    interaction_type=InteractionType.MJ,
    uniform_lambdas=[0.1, 0.5, 1.0, 2.0],
    vqe_max_iter=50,
)
analysis.run()   # uruchamia SamplingVQE dla każdego scenariusza
print(analysis.summary())
analysis.plot(output_dir=Path("output/field_analysis"))
```

Generuje 3 wykresy: `energy_vs_lambda.png`, `energy_comparison_bar.png`, `probability_distributions.png`.

**Znany problem:** etykiety z `λ` (U+03BB) crashują logger na Windows CP1250 — obliczenia się kończą, crash jest tylko w logu. Fix: zastąpić `λ=` przez `lam=` w etykietach loggera (nie w tytułach wykresów).

**35 testów** w `tests/test_field_influence_analysis.py`.

---

### Etap B1 — LigandBead ✅

**Plik:** [`src/particle/ligand_bead.py`](file:///c:\repos\quantum-protein-folding\src\particle\ligand_bead.py)

Ligand jako wolno poruszający się bead na kratce — bez turn-qubitów, wyłącznie qubity pozycji.

```python
from particle.ligand_bead import LigandBead, PositionEncoding

# BINARY (domyślny): ceil(log2(N)) qubitów
lig = LigandBead(symbol="L", index=0, num_lattice_nodes=8)
lig.num_position_qubits   # → 3

# Projektor na węzeł k (SparsePauliOp)
P_k = lig.position_projector(node_index=3)

# UNARY: N qubitów (one-hot)
lig2 = LigandBead("L", 0, num_lattice_nodes=4, encoding=PositionEncoding.UNARY)
lig2.num_position_qubits  # → 4
```

**Decyzja projektowa:** niezależna klasa (nie dziedziczy z `Bead`) — `Bead` wymaga turn-qubitów i `parent_chain_len`, ligand nie ma żadnego z tych elementów.

**49 testów** w `tests/test_ligand_bead.py`.

---

### Etap B2 — LigandInteraction ✅

**Plik:** [`src/interaction/ligand_interaction.py`](file:///c:\repos\quantum-protein-folding\src\interaction\ligand_interaction.py)

Model oddziaływania liganda z residuami aminokwasowymi. Dwa tryby:

```python
from interaction.ligand_interaction import LigandInteraction

# HP_LIKE: ligand jako H lub P (używa macierzy HP)
lig_hp = LigandInteraction.hp_like(ligand_hp_type="H")
lig_hp.get_energy("A")   # A jest hydrofobowy → -1.0
lig_hp.get_energy("R")   # R jest polarny → 0.0

# CUSTOM: własny słownik energii
lig_cu = LigandInteraction.custom(
    energy_map={"A": -2.5, "K": -0.3, "G": -0.1},
    default_energy=0.0,
)
lig_cu.get_energy("A")   # → -2.5
lig_cu.get_energy("W")   # → 0.0  (default, W nie w mapie)
```

**82 testy** w `tests/test_ligand_interaction.py`.

---

### Etap B3 — HamiltonianBuilder + H_ligand ✅

**Plik:** [`src/builder/hamiltonian_builder.py`](file:///c:\repos\quantum-protein-folding\src\builder\hamiltonian_builder.py)

Rozszerzona metoda `sum_hamiltonians()` o parametry liganda:

```python
H = builder.sum_hamiltonians(
    ligand=ligand_bead,                  # LigandBead
    ligand_interaction=ligand_interact,  # LigandInteraction
)
# → SparsePauliOp o rozmiarze (protein_qubits + n_pos) qubitów
```

**Architektura rejestru:**
```
[ turn-qubits (backbone) | ligand position qubits ]
   (N-1)*2 qubitów           ceil(log2(nodes)) qubitów
```

**Ważne — Wariant A:**  
`H_ligand = (Σ_i E_i) · I_total` — skalarne przesunięcie w pełnej przestrzeni. Ligand NIE jest przestrzennie sprzężony z białkiem. Rejestr qubitów liganda jest poprawnie dodany do układu (VQE widzi je), ale Hamiltonian nie narzuca preferowanej pozycji liganda.

**23 testy** w `tests/test_hamiltonian_ligand.py`.

---

### Etap B4 — LigandAnalysis ✅

**Plik:** [`src/analysis/ligand_analysis.py`](file:///c:\repos\quantum-protein-folding\src\analysis\ligand_analysis.py)

Pełna klasa analizy protein-ligand z VQE i wizualizacją.

```python
from analysis.ligand_analysis import LigandAnalysis
from enums import InteractionType

analysis = LigandAnalysis(
    main_chain="HPPHH",
    interaction_type=InteractionType.HP,
    ligand_hp_type="H",          # tryb HP_LIKE
    num_lattice_nodes=8,
    vqe_max_iter=100,
)
analysis.run()   # baseline + protein+ligand

# Energia wiązania: ΔE = E(protein+lig) - E(baseline)
delta_e = analysis.compute_binding_energy()

# Rozkład prawdopodobieństwa pozycji liganda
dist = analysis.compute_ligand_position_distribution()  # list[float], len=num_lattice_nodes

# Heurystyka enkapsulacji
enc = analysis.detect_encapsulation(interior_fraction=0.5, encapsulation_threshold=0.6)
print(enc.is_encapsulated, enc.score)

# Wykresy
analysis.plot(output_dir=Path("output/ligand"))
print(analysis.summary())
```

Dostępne gotowe funkcje eksperymentalne:
```python
from analysis.ligand_analysis import run_hp_hydrophobic_experiment, run_mj_strong_ligand_experiment

# Szybki eksperyment HP
hp = run_hp_hydrophobic_experiment("HPPHH", num_lattice_nodes=8, output_dir=Path("output/hp"))

# Eksperyment MJ z silnym ligandem
mj = run_mj_strong_ligand_experiment("ACDEFGH", num_lattice_nodes=8, output_dir=Path("output/mj"))
```

**55 testów** w `tests/test_ligand_analysis.py`.

---

## 4. Podsumowanie testów — aktualny stan

```
uv run pytest tests/ -q
312 passed, 22 warnings in ~19s
```

| Plik testów | Moduł | Testów | Status |
|-------------|-------|--------|--------|
| `test_external_field.py` | ExternalField (A1) | 41 | ✅ |
| `test_hamiltonian_external_field.py` | H_field (A2) | 13 | ✅ |
| `test_field_influence_analysis.py` | FieldInfluenceAnalysis (A3) | 35 | ✅ |
| `test_ligand_bead.py` | LigandBead (B1) | 49 | ✅ |
| `test_ligand_interaction.py` | LigandInteraction (B2) | 82 | ✅ |
| `test_hamiltonian_ligand.py` | H_ligand (B3) | 23 | ✅ |
| `test_ligand_analysis.py` | LigandAnalysis (B4) | 55 | ✅ |
| `test_utils.py` | Utilities | 14 | ✅ |
| **Łącznie** | | **312** | **✅** |

---

## 5. Znane ograniczenia (Wariant A)

> [!IMPORTANT]
> Wszystkie zaimplementowane człony Hamiltonianu (H_field i H_ligand) są **Wariantem A** — skalarnymi przesunięciami energetycznymi. Oznacza to:
>
> 1. **Optymalna konformacja białka NIE zmienia się** przy różnych wartościach λ (pole) lub różnych ligandach — bo różnice energii między stanami są zachowane.
> 2. **Ligand NIE jest przestrzennie sprzężony** z białkiem — jego qubity pozycji są w rejestrze, ale Hamiltonian nie preferuje żadnej konkretnej pozycji.
> 3. **Rozkład pozycji liganda** wynika z preferencji ansatzu VQE i parametrów CVaR, nie z fizycznego sprzężenia.
> 4. **Heurystyka enkapsulacji** jest 1D (indeks węzła jako proxy dla geometrii), nie 3D.

---

## 6. Co jeszcze musi zostać zrobione — Wariant B (priorytet!)

### Wariant B — Pełne przestrzenne sprzężenie liganda z białkiem

Wariant B to fundamentalna rozbudowa modelu. Wymaga przypisania jawnych qubitów pozycji KAŻDEMU beadowi białkowego łańcucha, co pozwoli na operator `δ(r_i = r_L)` — sprzężenie przestrzenne liganda z białkiem.

#### B-ext1: Qubity pozycji per bead białka

Aktualnie: beady białka mają tylko turn-qubity. Pozycja 3D wynika impliciite z sekwencji skrętów (oblicza ją `DistanceMap`).

Potrzeba: każdy bead `i` otrzymuje dodatkowy rejestr qubitów pozycji `P_bead_i^(k)` → projektor na węzeł `k`.

**Plik do stworzenia:** `src/particle/chain_bead_position.py`  
**Modele do zmodyfikowania:** `Bead`, `Protein`, `DistanceMap`

#### B-ext2: Operator δ(r_i = r_L)

Sprzężenie: jeśli bead `i` jest na węźle `k` ORAZ ligand jest na węźle `k`, dodaj energię kontaktu.

```
H_coupling = Σ_i Σ_k E_contact(aa_i) · P_bead_i^(k) ⊗ P_ligand^(k)
```

**Plik do zmodyfikowania:** `src/builder/hamiltonian_builder.py`  
Nowa metoda: `_build_ligand_contact_term_variant_b()`

#### B-ext3: Constraint — wykluczenie overlap protein–ligand

Ligand i bead białka nie mogą zajmować tego samego węzła jednocześnie.

```
H_overlap_penalty = Σ_k LARGE_CONST · P_bead_i^(k) ⊗ P_ligand^(k)  (suma po i)
```

#### B-ext4: Aktualizacja LigandAnalysis

Klasa `LigandAnalysis` musi zostać zaktualizowana aby:
- budować rozszerzony Hamiltonian (Wariant B)
- interpretować wyniki z większym rejestrem qubitów
- obsługiwać realną geometrię 3D przy heurystyce enkapsulacji

---

### Drobne poprawki do wykonania

#### Fix 1: Unicode w logerze FieldInfluenceAnalysis

**Plik:** `src/analysis/field_influence_analysis.py`  
**Problem:** etykiety `f"λ={lam} (uniform)"` crashują logger na Windows CP1250.  
**Fix:** zastąpić `λ=` przez `lam=` w wywołaniach `logger.*()` (nie w tytułach wykresów matplotlib — matplotlib używa własnych fontów).

Dotyczy linii z `logger.info(...)` gdzie pojawia się `λ`.

#### Fix 2: Deterministyczne wyniki VQE (seed)

Aktualnie VQE nie ma stałego ziarna — wyniki różnią się między uruchomieniami.  
**Fix:** dodać `seed_simulator=42` do `AerSimulator` w `src/backend/__init__.py` lub przekazać `seed` przez parametr.

#### Fix 3: Komentarze w pliku ligand_bead.py — PositionEncoding.BINARY vs UNARY

W metodzie `_build_position_qubits()` dla UNARY komentarz jest niespójny (linia 224 mówi o `half*(I - Z_k)` podczas gdy UNARY i BINARY używają tej samej funkcji `build_turn_qubit`). Warto dodać wyjaśniający komentarz.

---

### Następny krok rekomendowany

Jeśli nie chcesz zaczynać od Wariantu B (duży nakład pracy), możesz zacząć od:

1. **Notebook demonstracyjny dla B1-B4** — aktualizacja `docs/field_extension_presentation.ipynb` o etapy B1-B4. Aktualnie notebook pokrywa tylko etapy A1-A3 z syntetycznymi wynikami. Można go zaktualizować aby uruchamiał `LigandAnalysis` z prawdziwym VQE.

2. **Fix Unicode** (15 minut roboty, duże UX improvement).

3. **Wariant B** — B-ext1 jest blokerem dla wszystkiego dalszego.

---

## 7. Jak uruchomić

### Środowisko

```powershell
# Aktywacja środowiska uv (Windows PowerShell)
cd c:\repos\quantum-protein-folding
uv sync

# Uruchomienie testów
uv run pytest tests/ -q

# Uruchomienie main.py (pełny workflow dla APRLRFY)
uv run python src/main.py

# Uruchomienie FieldInfluenceAnalysis
uv run python -c "
import sys; sys.path.insert(0, 'src')
from analysis.field_influence_analysis import FieldInfluenceAnalysis
from enums import InteractionType
from pathlib import Path

a = FieldInfluenceAnalysis('HPPH', InteractionType.HP, uniform_lambdas=[0.1, 1.0], vqe_max_iter=30)
a.run()
print(a.summary())
"

# Uruchomienie LigandAnalysis
uv run python -c "
import sys; sys.path.insert(0, 'src')
from analysis.ligand_analysis import LigandAnalysis
from enums import InteractionType
from pathlib import Path

a = LigandAnalysis('HPPH', interaction_type=InteractionType.HP, ligand_hp_type='H', num_lattice_nodes=4, vqe_max_iter=30)
a.run()
print(a.summary())
print('Binding energy:', a.compute_binding_energy())
print('Position dist:', a.compute_ligand_position_distribution())
enc = a.detect_encapsulation()
print('Encapsulated:', enc.is_encapsulated, 'score:', round(enc.score, 4))
"
```

### Zależności (pyproject.toml)

```
qiskit >= 1.0
qiskit-nature
qiskit-algorithms
qiskit-aer
scipy
numpy
matplotlib
python-docx   # dodane w tej sesji
```

---

## 8. Kluczowe konwencje i decyzje projektowe

| Konwencja | Opis |
|-----------|------|
| **Import path** | W testach i skryptach: `sys.path.insert(0, 'src')`, potem `from particle.external_field import ...` (bez `src.`) |
| **Backward compat** | Wszystkie nowe parametry domyślnie `None` — `sum_hamiltonians()` bez argumentów działa identycznie jak przed |
| **TDD** | Każdy nowy moduł ma plik testów napisany równolegle z implementacją |
| **Styl** | Ruff + pyrefly type check (zob. `pyproject.toml`) |
| **Logging** | `from logger import get_logger; logger = get_logger()` — unikać `print()` w kodzie produkcyjnym |
| **SparsePauliOp** | Wszystkie operatory Hamiltonianu są typu `qiskit.quantum_info.SparsePauliOp` |
| **Wariant A** | Bieżąca implementacja — skalarne przesunięcia. Dokumentacja explicite mówi o tym ograniczeniu |

---

## 9. Pliki do przeczytania przed startem

W kolejności ważności:

1. [`walkthrough.md`](file:///c:\repos\quantum-protein-folding\walkthrough.md) — szczegółowy opis B1-B4 z diagramami
2. [`src/builder/hamiltonian_builder.py`](file:///c:\repos\quantum-protein-folding\src\builder\hamiltonian_builder.py) — docstring modułu wyjaśnia Wariant A matematycznie
3. [`src/analysis/ligand_analysis.py`](file:///c:\repos\quantum-protein-folding\src\analysis\ligand_analysis.py) — docstring modułu wyjaśnia architekturę pipeline
4. [`src/particle/ligand_bead.py`](file:///c:\repos\quantum-protein-folding\src\particle\ligand_bead.py) — design decision: dlaczego nie dziedziczy z Bead
5. [`src/main.py`](file:///c:\repos\quantum-protein-folding\src\main.py) — jak wygląda standardowy workflow end-to-end

---

## 10. Quick-start checklist dla nowego agenta

- [ ] Uruchom `uv run pytest tests/ -q` — powinno być `312 passed`
- [ ] Przeczytaj `walkthrough.md` — pełny opis etapów B1-B4
- [ ] Zdecyduj co robisz: Fix Unicode / Notebook / Wariant B
- [ ] Dla każdego nowego modułu: najpierw testy, potem implementacja (TDD)
- [ ] Sprawdzaj backward compat: `sum_hamiltonians()` bez argumentów musi dawać te same wyniki
- [ ] Po każdej zmianie: `uv run pytest tests/ -q` — wszystkie 312 muszą przejść
