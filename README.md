# Enthalpie

Calcule l'enthalpie (ΔH°), l'entropie (ΔS°), l'enthalpie libre (ΔG = ΔH − TΔS),
la constante d'équilibre (K = exp(−ΔG/RT)), la température d'inversion et la
spontanéité d'une réaction chimique à partir de données standard (298,15 K).

Aucune dépendance externe (Python 3.8+).

## Utilisation

```bash
python Enthalpie.py "2 H2 + O2 -> 2 H2O(g)"
python Enthalpie.py "CaCO3 -> CaO + CO2" -T 298.15 1200
python Enthalpie.py "N2 + 3 H2 -> 2 NH3" --plage 300 900 100
python Enthalpie.py --liste                       # espèces connues
python Enthalpie.py "A -> B" --donnees mes_donnees.json
```

- L'état physique peut être omis s'il est non ambigu (`H2O` demande `(l)` ou `(g)`).
- L'équilibrage de l'équation est vérifié (avertissement sinon).
- Fichier `--donnees` : `{"X(g)": [ΔfH en kJ/mol, S en J/mol/K]}`.

## En Python

```python
from Enthalpie import ChemicalReaction
r = ChemicalReaction.from_equation("CaCO3 -> CaO + CO2")
r.calculate_delta_H(); r.calculate_delta_G(1200); r.inversion_temperature()
```

## Hypothèse

ΔH° et ΔS° sont supposés indépendants de la température (approximation d'Ellingham) ;
les résultats loin de 298 K sont donc approximatifs.

## Tests

```bash
python -m unittest
```
