class ChemicalReaction:
    def __init__(self, reactants, products, enthalpies, entropies):
        self.reactants = reactants  # Dictionnaire des réactifs et leurs coefficients stœchiométriques
        self.products = products    # Dictionnaire des produits et leurs coefficients stœchiométriques
        self.enthalpies = enthalpies  # Dictionnaire des enthalpies de formation (kJ/mol)
        self.entropies = entropies    # Dictionnaire des entropies (J/mol*K)

    def calculate_delta_H(self):
        delta_H = sum(self.products[compound] * self.enthalpies[compound] for compound in self.products) - \
                  sum(self.reactants[compound] * self.enthalpies[compound] for compound in self.reactants)
        return delta_H

    def calculate_delta_S(self):
        delta_S = sum(self.products[compound] * self.entropies[compound] for compound in self.products) - \
                  sum(self.reactants[compound] * self.entropies[compound] for compound in self.reactants)
        return delta_S

    def calculate_delta_G(self, temperature):
        delta_H = self.calculate_delta_H()
        delta_S = self.calculate_delta_S()
        delta_G = delta_H - (temperature * delta_S / 1000)  # Convertir J en kJ pour delta_S
        return delta_G

def main():
    # Exemple de réaction : 2 H2 + O2 -> 2 H2O
    reactants = {'H2': 2, 'O2': 1}
    products = {'H2O': 2}
    enthalpies = {'H2': 0, 'O2': 0, 'H2O': -241.8}  # Enthalpies de formation standard (kJ/mol)
    entropies = {'H2': 130.6, 'O2': 205.2, 'H2O': 188.7}  # Entropies standard (J/mol*K)

    reaction = ChemicalReaction(reactants, products, enthalpies, entropies)

    delta_H = reaction.calculate_delta_H()
    delta_S = reaction.calculate_delta_S()
    temperature = 298.15  # Température standard en Kelvin
    delta_G = reaction.calculate_delta_G(temperature)

    print(f"Delta H (Enthalpie) : {delta_H:.2f} kJ/mol")
    print(f"Delta S (Entropie) : {delta_S:.2f} J/mol*K")
    print(f"Delta G (Enthalpie libre) à {temperature} K : {delta_G:.2f} kJ/mol")

if __name__ == "__main__":
    main()