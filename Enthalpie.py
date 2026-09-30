#!/usr/bin/env python3
"""Thermodynamique d'une réaction chimique : ΔH, ΔS, ΔG, K et spontanéité.

Exemples :
    python Enthalpie.py "2 H2 + O2 -> 2 H2O(g)"
    python Enthalpie.py "CaCO3(s) -> CaO(s) + CO2(g)" -T 298.15 1000 1200
    python Enthalpie.py "N2 + 3 H2 -> 2 NH3" --plage 300 900 100
    python Enthalpie.py --liste

Hypothèse : ΔH° et ΔS° sont indépendants de la température (approximation
d'Ellingham), valeurs tabulées à 298,15 K.
"""
import argparse
import json
import math
import re
import sys

R = 8.314462618  # Constante des gaz parfaits (J/mol/K)
T_STANDARD = 298.15  # K

# Données standard à 298,15 K : espèce -> (ΔfH° en kJ/mol, S° en J/mol/K)
COMPOUNDS = {
    "H2(g)": (0.0, 130.68),
    "O2(g)": (0.0, 205.15),
    "N2(g)": (0.0, 191.61),
    "Cl2(g)": (0.0, 223.07),
    "C(graphite)": (0.0, 5.74),
    "Fe(s)": (0.0, 27.28),
    "Al(s)": (0.0, 28.33),
    "H2O(l)": (-285.83, 69.95),
    "H2O(g)": (-241.83, 188.84),
    "CO(g)": (-110.53, 197.66),
    "CO2(g)": (-393.51, 213.79),
    "CH4(g)": (-74.6, 186.26),
    "C2H6(g)": (-84.0, 229.2),
    "C2H4(g)": (52.4, 219.3),
    "C3H8(g)": (-103.85, 270.3),
    "C6H6(l)": (49.1, 173.4),
    "CH3OH(l)": (-239.2, 126.8),
    "C2H5OH(l)": (-277.6, 160.7),
    "C6H12O6(s)": (-1273.3, 212.1),
    "NH3(g)": (-45.9, 192.77),
    "NO(g)": (90.25, 210.76),
    "NO2(g)": (33.18, 240.06),
    "N2O4(g)": (9.16, 304.29),
    "HCl(g)": (-92.31, 186.9),
    "H2S(g)": (-20.6, 205.8),
    "SO2(g)": (-296.83, 248.22),
    "SO3(g)": (-395.72, 256.76),
    "CaCO3(s)": (-1206.9, 92.9),
    "CaO(s)": (-635.09, 38.1),
    "NaCl(s)": (-411.15, 72.13),
    "Fe2O3(s)": (-824.2, 87.4),
    "Al2O3(s)": (-1675.7, 50.92),
}


def _strip_state(species):
    return re.sub(r"\((?:g|l|s|aq|graphite)\)$", "", species)


def element_counts(species):
    """Compte les atomes d'une formule ('Fe2O3(s)', 'Ca(OH)2'...)."""
    formula = _strip_state(species)
    if species == "C(graphite)":
        formula = "C"
    stack = [{}]
    pos = 0
    token = re.compile(r"([A-Z][a-z]?)(\d*)|(\()|\)(\d*)")
    while pos < len(formula):
        m = token.match(formula, pos)
        if not m:
            raise ValueError(f"Formule invalide : '{species}'")
        pos = m.end()
        if m.group(1):
            n = int(m.group(2) or 1)
            stack[-1][m.group(1)] = stack[-1].get(m.group(1), 0) + n
        elif m.group(3):
            stack.append({})
        else:
            if len(stack) == 1:
                raise ValueError(f"Parenthèses déséquilibrées : '{species}'")
            n = int(m.group(4) or 1)
            inner = stack.pop()
            for e, c in inner.items():
                stack[-1][e] = stack[-1].get(e, 0) + c * n
    if len(stack) != 1:
        raise ValueError(f"Parenthèses déséquilibrées : '{species}'")
    return stack[0]


def resolve_species(name, database):
    """Retourne la clé de la base pour 'name' (l'état physique peut être omis)."""
    if name in database:
        return name
    candidates = [k for k in database if _strip_state(k) == name or
                  (k.startswith(name + "(") and k.endswith(")"))]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise KeyError(f"Espèce inconnue : '{name}' (voir --liste ou --donnees)")
    raise KeyError(f"'{name}' est ambigu, précisez l'état : {', '.join(candidates)}")


def _parse_side(side):
    result = {}
    for term in side.split("+"):
        term = term.strip()
        if not term:
            raise ValueError("Terme vide dans l'équation")
        m = re.fullmatch(r"(\d+(?:[.,]\d+)?)?\s*(\S.*)", term)
        coeff = float(m.group(1).replace(",", ".")) if m.group(1) else 1.0
        if coeff <= 0:
            raise ValueError(f"Coefficient invalide dans '{term}'")
        coeff = int(coeff) if coeff == int(coeff) else coeff
        name = m.group(2).strip()
        result[name] = result.get(name, 0) + coeff
    return result


def parse_equation(equation):
    """'2 H2 + O2 -> 2 H2O' -> ({'H2': 2, 'O2': 1}, {'H2O': 2})."""
    parts = re.split(r"->|→|=", equation)
    if len(parts) != 2:
        raise ValueError("L'équation doit contenir une seule flèche '->'")
    return _parse_side(parts[0]), _parse_side(parts[1])


class ChemicalReaction:
    def __init__(self, reactants, products, enthalpies, entropies):
        """reactants/products : {espèce: coefficient stœchiométrique}
        enthalpies : {espèce: ΔfH° en kJ/mol} ; entropies : {espèce: S° en J/mol/K}
        """
        if not reactants or not products:
            raise ValueError("Une réaction nécessite des réactifs et des produits")
        for species in list(reactants) + list(products):
            for label, table in (("enthalpie", enthalpies), ("entropie", entropies)):
                if species not in table:
                    raise KeyError(f"{label} manquante pour '{species}'")
        self.reactants = reactants
        self.products = products
        self.enthalpies = enthalpies
        self.entropies = entropies

    @classmethod
    def from_equation(cls, equation, database=None):
        database = COMPOUNDS if database is None else database
        raw_r, raw_p = parse_equation(equation)
        reactants = {resolve_species(k, database): v for k, v in raw_r.items()}
        products = {resolve_species(k, database): v for k, v in raw_p.items()}
        enthalpies = {k: v[0] for k, v in database.items()}
        entropies = {k: v[1] for k, v in database.items()}
        return cls(reactants, products, enthalpies, entropies)

    def _delta(self, table):
        return (sum(n * table[c] for c, n in self.products.items())
                - sum(n * table[c] for c, n in self.reactants.items()))

    def imbalance(self):
        """Atomes en excès (produits - réactifs) ; {} si l'équation est équilibrée."""
        totals = {}
        for side, sign in ((self.products, 1), (self.reactants, -1)):
            for species, n in side.items():
                for elem, c in element_counts(species).items():
                    totals[elem] = totals.get(elem, 0) + sign * n * c
        return {e: v for e, v in totals.items() if abs(v) > 1e-9}

    def is_balanced(self):
        return not self.imbalance()

    def calculate_delta_H(self):
        """ΔH° de réaction (kJ/mol)."""
        return self._delta(self.enthalpies)

    def calculate_delta_S(self):
        """ΔS° de réaction (J/mol/K)."""
        return self._delta(self.entropies)

    def calculate_delta_G(self, temperature=T_STANDARD):
        """ΔG = ΔH - TΔS (kJ/mol), température en kelvin."""
        if temperature <= 0:
            raise ValueError("La température doit être > 0 K")
        return self.calculate_delta_H() - temperature * self.calculate_delta_S() / 1000

    def equilibrium_constant(self, temperature=T_STANDARD):
        """K = exp(-ΔG / RT) ; float('inf') en cas de dépassement."""
        exponent = -self.calculate_delta_G(temperature) * 1000 / (R * temperature)
        try:
            return math.exp(exponent)
        except OverflowError:
            return float("inf")

    def inversion_temperature(self):
        """T (K) où ΔG = 0, ou None si ΔH et ΔS sont de même signe nul/opposé
        (pas de changement de spontanéité pour T > 0)."""
        dH, dS = self.calculate_delta_H(), self.calculate_delta_S()
        if dS == 0:
            return None
        t = dH * 1000 / dS
        return t if t > 0 else None

    def spontaneity(self, temperature=T_STANDARD):
        dG = self.calculate_delta_G(temperature)
        if abs(dG) < 1e-9:
            return "à l'équilibre"
        return "spontanée" if dG < 0 else "non spontanée"

    def describe_type(self):
        dH = self.calculate_delta_H()
        return "exothermique" if dH < 0 else "endothermique" if dH > 0 else "athermique"

    def equation(self):
        def side(d):
            return " + ".join(f"{n if n != 1 else ''}{' ' if n != 1 else ''}{s}"
                              for s, n in d.items())
        return f"{side(self.reactants)} -> {side(self.products)}"


def format_K(k):
    return "∞" if math.isinf(k) else f"{k:.3e}"


def report(reaction, temperatures):
    lines = [f"Réaction : {reaction.equation()}"]
    imb = reaction.imbalance()
    if imb:
        lines.append("ATTENTION - équation non équilibrée : " +
                     ", ".join(f"{e} ({v:+g})" for e, v in imb.items()))
    dH, dS = reaction.calculate_delta_H(), reaction.calculate_delta_S()
    lines.append(f"Delta H° (Enthalpie)  : {dH:.2f} kJ/mol ({reaction.describe_type()})")
    lines.append(f"Delta S° (Entropie)   : {dS:.2f} J/mol/K")
    t_inv = reaction.inversion_temperature()
    if t_inv:
        lines.append(f"Température d'inversion (ΔG = 0) : {t_inv:.1f} K")
    for t in temperatures:
        dG = reaction.calculate_delta_G(t)
        lines.append(f"Delta G à {t:g} K (Enthalpie libre) : {dG:.2f} kJ/mol "
                     f"-> {reaction.spontaneity(t)}, K = {format_K(reaction.equilibrium_constant(t))}")
    return "\n".join(lines)


def load_database(path):
    """Fichier JSON : {"espèce": [ΔfH kJ/mol, S J/mol/K], ...}"""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    db = dict(COMPOUNDS)
    for k, v in data.items():
        if not (isinstance(v, (list, tuple)) and len(v) == 2):
            raise ValueError(f"Entrée invalide pour '{k}' : attendu [ΔfH, S]")
        db[k] = (float(v[0]), float(v[1]))
    return db


def build_parser():
    p = argparse.ArgumentParser(
        description="Calcule ΔH, ΔS, ΔG, K et la spontanéité d'une réaction.")
    p.add_argument("equation", nargs="?", help='ex. "2 H2 + O2 -> 2 H2O(g)"')
    p.add_argument("-T", "--temperature", type=float, nargs="+", default=None,
                   help="température(s) en K (défaut : 298.15)")
    p.add_argument("--plage", type=float, nargs=3, metavar=("TMIN", "TMAX", "PAS"),
                   help="balayage de températures")
    p.add_argument("--donnees", metavar="FICHIER.json",
                   help="données thermodynamiques supplémentaires")
    p.add_argument("--liste", action="store_true", help="liste les espèces connues")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        database = load_database(args.donnees) if args.donnees else COMPOUNDS
        if args.liste:
            print(f"{'Espèce':<14}{'ΔfH° (kJ/mol)':>15}{'S° (J/mol/K)':>15}")
            for k, (h, s) in sorted(database.items()):
                print(f"{k:<14}{h:>15.2f}{s:>15.2f}")
            return 0
        equation = args.equation or "2 H2 + O2 -> 2 H2O(g)"
        temps = list(args.temperature or [])
        if args.plage:
            tmin, tmax, pas = args.plage
            if pas <= 0 or tmax < tmin:
                raise ValueError("--plage : PAS > 0 et TMAX >= TMIN requis")
            t = tmin
            while t <= tmax + 1e-9:
                temps.append(t)
                t += pas
        if not temps:
            temps = [T_STANDARD]
        print(report(ChemicalReaction.from_equation(equation, database), temps))
        return 0
    except (ValueError, KeyError, OSError) as exc:
        print(f"Erreur : {exc.args[0] if isinstance(exc, KeyError) else exc}",
              file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
