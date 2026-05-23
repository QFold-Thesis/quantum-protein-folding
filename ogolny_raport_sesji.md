# Raport ogólny — Rozszerzenie modelu kwantowego zwijania białka

> **Projekt:** Kwantowe symulacje struktury białka na kratce FCC/diamentowej  
> **Sesja:** Wariant A — Zewnętrzne pole oddziaływań  
> **Etapy ukończone:** A1 · A2 · A3

---

## Co zostało zbudowane

W tej sesji rozszerzono istniejący model kwantowy zwijania białka o możliwość uwzględnienia **zewnętrznego pola oddziaływań** — dodatkowego składnika energetycznego, który modeluje obecność ligandu, środowiska membranowego lub aktywnego centrum katalitycznego.

Implementacja obejmuje trzy warstwy:

```
ExternalField  ──▶  HamiltonianBuilder  ──▶  FieldInfluenceAnalysis
(definicja pola)    (kwantowy Hamiltonian)   (analiza porównawcza VQE)
```

---

## Interpretacja fizyczna

### Model bazowy

Białko jest reprezentowane jako łańcuch beadów na kratce FCC (tetrahedralnej). Hamiltoniam kwantowy:

$$H = H_{\text{backbone}} + H_{\text{backtrack}}$$

opisuje energię wewnętrzną łańcucha głównego i penalizuje niedozwolone skręty. VQE (Variational Quantum Eigensolver) minimalizuje energię i znajduje optymalną konformację.

### Rozszerzenie — Wariant A

Dodano człon pola zewnętrznego:

$$H = H_{\text{backbone}} + H_{\text{backtrack}} + H_{\text{field}}$$

gdzie:

$$H_{\text{field}} = \sum_{i=0}^{N-1} E_{\text{field}}(i) \cdot \mathbf{I}$$

Każdy bead $i$ wnosi do Hamiltonianu energię zależną od jego pozycji w sekwencji. Pole może być:
- **jednorodne** — stała wartość $\lambda$ dla każdego beadu (testowanie globalnej stabilizacji)
- **niejednorodne** — różne wartości dla różnych pozycji (np. silniejsze pole w centrum, modelujące kieszeń wiążącą)

### Ograniczenia i interpretacja

Wariant A jest **polem sekwencyjnym** — przypisuje energię beadom według ich numeru w łańcuchu, nie według faktycznej pozycji przestrzennej. Oznacza to, że pole jednorodne przesuwa energię wszystkich konformacji o stałą wartość, a pole niejednorodne różnicuje je według sekwencji.

**Pełne sprzężenie przestrzenne** (tj. pole o sile zależnej od tego, gdzie faktycznie bead *wyląduje* na kratce po optymalizacji kwantowej) wymaga Wariantu B — liganda jako osobnego beadu z własnymi stopniami swobody kwantowymi.

---

## Co umożliwia nowy kod

### 1. Definiowanie pola

```python
from particle.external_field import ExternalField

# Jednorodne pole o sile -1.0
pole = ExternalField.uniform(strength=-1.0)

# Niejednorodne — wzmocnienie w centrum łańcucha 7-beadowego
pole = ExternalField.non_uniform({
    (0,): -0.1,  (1,): -0.5,  (2,): -0.9,
    (3,): -1.0,  (4,): -0.9,  (5,): -0.5,  (6,): -0.1,
})
```

### 2. Budowanie Hamiltonianu z polem

```python
from builder.hamiltonian_builder import HamiltonianBuilder

builder = HamiltonianBuilder(
    protein=protein,
    interaction=interaction,
    distance_map=distance_map,
    contact_map=contact_map,
    external_field=pole,   # ← nowy parametr; None = jak poprzednio
)
H = builder.sum_hamiltonians()
```

### 3. Analiza porównawcza

```python
from analysis.field_influence_analysis import FieldInfluenceAnalysis

analiza = FieldInfluenceAnalysis(
    main_chain="HPPHH",
    uniform_lambdas=[0.1, 0.5, 1.0, 2.0],
    vqe_max_iter=100,
)
analiza.run()
print(analiza.summary())
analiza.plot(output_dir=Path("wyniki/"))
```

Automatycznie generowane wykresy:
- **Energia minimum vs. λ** — jak silne pole zmienia najgłębsze minimum energetyczne
- **Porównanie słupkowe** — energia we wszystkich scenariuszach
- **Rozkłady prawdopodobieństwa** — top-8 stanów per scenariusz

---

## Wyniki ilościowe

| Metryka | Wartość |
|---------|---------|
| Nowe pliki źródłowe | 4 |
| Zmodyfikowane pliki | 2 |
| Nowe testy (szybkie) | 101 |
| Nowe testy (slow/e2e) | 1 |
| Łączny czas testów | < 8 s |
| Regresje w istniejących testach | 0 |

---

## Co dalej — Wariant B

Wariant B zakłada modelowanie ligandu jako **pełnoprawnego beadu** z własną pozycją kwantową:

- **B1:** `LigandBead` — klasa beadu z pozycją zakodowaną w qubitach
- **B2:** Operator `H_ligand-chain` = prawdziwy operator kontaktu $\delta(r_i, r_{\text{lig}})$
- **B3:** Wspólna optymalizacja VQE dla białko + ligand
- **B4:** Analiza porównawcza krajobrazu energetycznego z ligandem vs. bez

---

*Etapy A1–A3 kompletne. Kod zweryfikowany testami, gotowy do dalszego rozszerzenia.*
