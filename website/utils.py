
def render_mol_formula(molformula):
    chars = [char for char in molformula]
    for char in chars:
        if char.isnumeric():
            chars[chars.index(char)] = "<sub>" + char + "</sub>"
    return "".join(chars)

