# Raport szczegółowy — Sesja implementacyjna
## Rozszerzenie modelu zwijania białka o zewnętrzne pole oddziaływań

> **Projekt:** QFold-Thesis / quantum-protein-folding  
> **Data sesji:** 2026-05-23  
> **Wariant:** A — Zewnętrzne pole oddziaływań (etapy A1 -> A2 -> A3)  
> **Status końcowy:** ✅ 102 testy przechodzą w < 8 s (bez wolnych testów end-to-end)

---

## 1. Kontekst i cel

Celem sesji było rozszerzenie istniejącego kwantowego modelu zwijania białka (representacja FCC/diamentowa, oddziaływania HP i MJ) o **zewnętrzne pole oddziaływań** (Wariant A). Implementacja podzielona była na trzy etapy:

| Etap | Nazwa | Cel |
|------|-------|-----|
| A1 | `ExternalField` | Enkapsulacja pola jako obiektu z mapą `coord -> energia` |
| A2 | Integracja z `HamiltonianBuilder` | Dodanie członu `H_field` do Hamiltonianu |
| A3 | `FieldInfluenceAnalysis` | Analiza porównawcza i wizualizacja wyników VQE |

---

## 2. Etap A1 — Klasa `ExternalField`

### 2.1 Lokalizacja

```
src/particle/__init__.py          ← nowy pakiet
src/particle/external_field.py    ← implementacja
tests/test_external_field.py      ← 41 testów jednostkowych
```

### 2.2 Architektura klasy

```python
class ExternalField:
    mode: FieldMode          # UNIFORM | NON_UNIFORM
    default_energy: float    # fallback dla brakujących węzłów
    _energy_map: dict[tuple[int, ...], float]

    # Fabryki
    @classmethod def uniform(strength: float) -> ExternalField
    @classmethod def non_uniform(energy_map, default_energy=0.0) -> ExternalField

    # API
    def get_energy(lattice_coords: tuple[int, ...]) -> float
    def set_energy(lattice_coords, energy: float) -> None
    def nodes() -> dict   # kopia mapy
```

### 2.3 Kluczowe decyzje projektowe

**Tryb UNIFORM:** `default_energy = strength`, `_energy_map = {}`. Wywołanie `get_energy()` korzysta z `_energy_map.get(coords, default_energy)` — umożliwia to nadpisywanie węzłów przez `set_energy()` nawet w trybie jednorodnym.

**Tryb NON_UNIFORM:** Dowolna mapa `(i, j, ...)->float`; brakujące węzły zwracają `default_energy` (domyślnie 0.0).

**Walidacja:** Każda metoda sprawdza:
- `lattice_coords` musi być `tuple` (nie `list`, nie `int`)
- Wszystkie wartości energii muszą być skończone (`math.isfinite`)
- `energy_map` musi być słownikiem

**Agnostyczność geometryczna:** Klucze mogą być krotek dowolnego wymiaru — 2D, 3D, tetrahedralne. Nie jest zakodowana żadna konkretna geometria kratki.

### 2.4 Naprawa błędu znalezionego w testach

Pierwotna implementacja `get_energy` w trybie UNIFORM ignorowała słownik i zwracała zawsze `default_energy`, co uniemożliwiało `set_energy()` na polach jednorodnych. Poprawka: ujednolicenie ścieżki — zawsze `_energy_map.get(coords, default_energy)`.

### 2.5 Wyniki testów (A1)

```
41 passed in 0.42s
Klasy testów:
  TestUniformFactory (9)         — fabryka, tryb, wartości
  TestNonUniformFactory (11)     — fabryka, mapa, walidacja
  TestGetEnergyValidation (6)    — TypeError, wymiary koordynat
  TestSetEnergy (6)              — nadpisywanie, walidacja
  TestNodes (3)                  — kopiowanie, niezależność
  TestRepr (3)                   — format reprezentacji
  TestDirectConstruction (3)     — bezpośredni __init__
```

---

## 3. Etap A2 — Integracja z `HamiltonianBuilder`

### 3.1 Lokalizacja zmian

```
src/builder/hamiltonian_builder.py   ← zmodyfikowany (nie nowy)
tests/test_hamiltonian_external_field.py  ← 13 testów integracyjnych
```

### 3.2 Zmiany w `HamiltonianBuilder`

#### `__init__` — nowy parametr

```python
def __init__(
    self,
    protein: Protein,
    interaction: Interaction,
    distance_map: DistanceMap,
    contact_map: ContactMap,
    external_field: ExternalField | None = None,   # ← nowy
) -> None:
```

Wartość domyślna `None` gwarantuje **pełną kompatybilność wsteczną** — żaden istniejący kod wywołujący konstruktor nie wymaga zmian.

#### `sum_hamiltonians` — rozszerzenie

```python
h_backbone  = self._build_backbone_contact_term()
h_backtrack = self._add_backtracking_penalty()
h_field     = self._build_external_field_term()   # ← nowe wywołanie

part_hamiltonians = [h_backbone, h_backtrack, h_field]
```

#### `_build_external_field_term()` — nowa metoda

```python
def _build_external_field_term(self) -> SparsePauliOp:
    n_turn_qubits = (len(self.protein.main_chain) - 1) * QUBITS_PER_TURN

    if self.external_field is None:
        return build_identity_op(n_turn_qubits, EMPTY_OP_COEFF)   # zero-op

    h_field = build_identity_op(n_turn_qubits, EMPTY_OP_COEFF)    # startowe zero

    for bead in self.protein.main_chain:
        bead_coord = (bead.index,)
        energy = self.external_field.get_energy(bead_coord)
        h_field = h_field + energy * build_identity_op(n_turn_qubits)

    return h_field.simplify()
```

### 3.3 Formalizm fizyczny (Wariant A)

$$H_{\text{field}} = \sum_{i=0}^{N-1} E_{\text{field}}((i,)) \cdot \mathbf{I}$$

Gdzie:
- $i$ — sekwencyjny indeks beadu (koordynata kratki)
- $E_{\text{field}}$ — energia pola w danym węźle (zwracana przez `ExternalField.get_energy`)
- $\mathbf{I}$ — operator tożsamości na rejestrze qubitów skrętów

**Dlaczego nie pełne `δ(r_i, r_field)`?**  
W tym modelu pozycje beadów na kratce FCC są wynikiem *klasycznej* interpretacji sekwencji skrętów kwantowych (sumowanie wektorów bazowych FCC). Nie istnieje operator kwantowy pozycji `r_i` — jest on funkcją wszystkich qubitów skrętów i obliczalny dopiero po pomiarze. Pełne sprzężenie przestrzenne wymaga reprezentacji pozycji jako kwantowego stopnia swobody — to jest zadanie Wariantu B (bead ligandowy).

### 3.4 Weryfikacja kompatybilności wstecznej

Test `test_sum_hamiltonians_identical_without_field` weryfikuje:

```python
diff = (H_with_none_field - H_baseline).simplify()
for coeff in diff.coeffs:
    assert abs(coeff) < 1e-10   # dokładnie zero
```

### 3.5 Wyniki testów (A2)

```
13 passed in 0.78s
Klasy testów:
  TestNoneFieldBackwardCompatibility (3) — none=baseline, qubit count
  TestBuildExternalFieldTerm (5)         — zero-op, skalowanie I, liczba qubitów
  TestSumHamiltonians (5)                — shift Δ=N·s·I, qubit count, storage
```

---

## 4. Etap A3 — `FieldInfluenceAnalysis`

### 4.1 Lokalizacja

```
src/analysis/__init__.py                   ← nowy pakiet
src/analysis/field_influence_analysis.py   ← implementacja (~430 linii)
tests/test_field_influence_analysis.py     ← 35 testów (34 szybkie + 1 slow)
```

### 4.2 Architektura klasy

```
FieldInfluenceAnalysis
  ├── __init__(main_chain, interaction_type, uniform_lambdas, vqe_max_iter, ...)
  │     └── inicjalizuje Protein, ContactMap, DistanceMap, Interaction (reużywane)
  │
  ├── run()
  │     ├── _run_scenario("baseline", field=None)
  │     ├── _run_scenario("λ=0.1 (uniform)", ExternalField.uniform(0.1))
  │     ├── ...                                                ← sweep λ
  │     └── _run_scenario("non-uniform (centre boost)", _build_non_uniform_field())
  │
  ├── plot(output_dir)
  │     ├── _plot_energy_vs_lambda()    -> energy_vs_lambda.png
  │     ├── _plot_energy_bar()          -> energy_comparison_bar.png
  │     └── _plot_probability_distributions() -> probability_distributions.png
  │
  └── summary() -> str (tabela tekstowa)
```

### 4.3 `ScenarioResult` — kontener danych

```python
@dataclass
class ScenarioResult:
    label: str
    field: ExternalField | None
    minimum_energy: float
    best_bitstring: str
    state_probabilities: dict[str, float]   # bitstring -> prawdopodobieństwo
    vqe_iterations: list[int]
    vqe_energies: list[float]
```

### 4.4 Profil pola niejednorodnego

Funkcja Gaussa wycentrowana na środku łańcucha:

$$E_{\text{field}}((i,)) = \lambda_{\max} \cdot \exp\!\left(-\frac{(i - \mu)^2}{\sigma^2}\right)$$

Gdzie $\mu = (N-1)/2$, $\sigma = N/4$, $\lambda_{\max} = -1.0$ (atrakcyjne).

Efekt: silne pole w centrum łańcucha (aktywne centrum katalityczne), zanikające do ~2% na końcach.

### 4.5 Wizualizacja — ciemny motyw

Wszystkie wykresy używają spójnego stylu:
- Tło figury: `#0f1117`, osie: `#1a1d27`
- Paleta Matplotlib Category10 — konsekwentna między wykresami
- Siatka major + minor, ramka `#333`
- Zapisane jako PNG 150 dpi z `bbox_inches='tight'`

### 4.6 Zarządzanie testami

Rejestracja markera `slow` w `pyproject.toml`:

```toml
[tool.pytest.ini_options]
pythonpath = ["src"]
markers = [
    "slow: marks tests as slow (VQE end-to-end); deselect with '-m \"not slow\"'",
]
```

Uruchamianie:
```bash
pytest tests/ -m "not slow"   # 102 testy, ~8s
pytest tests/ -m slow         # pełny end-to-end (kilka minut)
```

### 4.7 Wyniki testów (A3)

```
34 passed (fast), 1 deselected (slow) in 6.70s
Klasy testów:
  TestInit (8)               — konfiguracja, łańcuchy, tryby
  TestBuildNonUniformField (5) — Gauss: max w centrum, skończone energię
  TestSummary (3)             — guard, zawartość, format
  TestPlotGuard (1)           — guard przed run()
  TestPlotOutput (5)          — PNG zapisane, katalog tworzony
  TestScenarioResult (4)      — dataclass, pola
  TestRunOrchestration (8)    — mock VQE: kolejność, etykiety, typy pól
  TestEndToEnd (1, slow)      — pełny pipeline
```

---

## 5. Naprawione błędy i problemy techniczne

| Problem | Diagnoza | Rozwiązanie |
|---------|----------|-------------|
| `get_energy` ignorował mapę w trybie UNIFORM | `if mode == UNIFORM: return default` — gałąź omijała słownik | Ujednolicenie: zawsze `_map.get(k, default)` |
| Testy przekraczały limit czasu przy `to_matrix()` | `2^n` macierz dla `n=20` qubitów zajmuje ~8 GB | Zamiana na porównanie współczynników Pauliego |
| `isinstance(ExternalField, ExternalField) == False` | Podwójny import przez `src.particle.*` vs `particle.*` — dwa obiekty klasy w cache | Usunięcie prefiksu `src.` z importów testowych |
| Błąd Tkinter `TclError` w testach wykresów | `matplotlib` próbuje otworzyć GUI bez serwera X | `os.environ['MPLBACKEND'] = 'Agg'` przed importem |

---

## 6. Całkowite statystyki testów po sesji

```
tests/test_external_field.py              41 passed
tests/test_hamiltonian_external_field.py  13 passed
tests/test_field_influence_analysis.py    34 passed (fast) + 1 slow
tests/test_utils.py                       14 passed (bez zmian, brak regresji)

RAZEM (fast):  102 passed, 1 deselected in 7.56s
```

---

## 7. Drzewo nowych i zmodyfikowanych plików

```
quantum-protein-folding/
├── src/
│   ├── particle/                          [NEW]
│   │   ├── __init__.py                    [NEW] ← eksportuje ExternalField
│   │   └── external_field.py              [NEW] ← 321 linii
│   ├── builder/
│   │   └── hamiltonian_builder.py         [MODIFIED] ← +external_field param
│   └── analysis/                          [NEW]
│       ├── __init__.py                    [NEW] ← eksportuje FieldInfluenceAnalysis
│       └── field_influence_analysis.py    [NEW] ← ~430 linii
├── tests/
│   ├── test_external_field.py             [NEW] ← 41 testów
│   ├── test_hamiltonian_external_field.py [NEW] ← 13 testów
│   └── test_field_influence_analysis.py   [NEW] ← 35 testów
└── pyproject.toml                         [MODIFIED] ← dodano marker 'slow'
```

---

## 8. Ograniczenia i przyszłe kierunki

### Ograniczenia Wariantu A

1. `H_field = Σ E_field((i,)) · I` jest **operatorem skalarnym** — nie rozkłada stanów kwantowych (jednakowo przesuwa energię wszystkich konformacji).
2. Nie ma **sprzężenia geometrycznego**: pole nie wie, gdzie faktycznie jest bead na kratce.
3. Pole niejednorodne różnicuje energię tylko **po sekwencji**, nie po pozycji przestrzennej.

### Następne kroki — Wariant B

- **B1:** `LigandBead` — nowy typ beadu z zakodowaną pozycją jako qubity
- **B2:** Operator `H_ligand` z prawdziwym `δ(r_i, r_ligand)` — kontakt bead↔ligand jako operator Pauliego
- **B3:** Optymalizacja pozycji ligandu i konformacji białka jednocześnie
- **B4:** Analiza porównawcza A vs B

---

*Wygenerowano automatycznie na podstawie historii sesji.*
