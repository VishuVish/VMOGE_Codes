import numpy as np

# Define Jones matrix elements as strings
Jxx, Jxy, Jyx, Jyy = 'J_{xx}', 'J_{xy}', 'J_{yx}', 'J_{yy}'
J_elements = [Jxx, Jxy, Jyx, Jyy]
J_conj_elements = [f'{elem}^*' for elem in J_elements]

# Kronecker product J ⊗ J* symbolically (as strings)
J_kron_Jconj = np.array([f"{j}*{jc}" for j in J_elements for jc in J_conj_elements]).reshape(4, 4)

# Transformation matrix A and its inverse
A = np.array([
    [1, 1, 0, 0],
    [1, -1, 0, 0],
    [0, 0, 1, 1],
    [0, 0, i, -i]
])
A_inv = np.linalg.inv(A)

# Compute Mueller matrix M = A (J ⊗ J*) A^(-1) symbolically
M = np.zeros((4, 4), dtype=object)
for i in range(4):
    for j in range(4):
        expr = []
        for k in range(4):
            for l in range(4):
                coeff = A[i, k] * A_inv[l, j]
                if abs(coeff) > 1e-10:  # Check for non-zero coefficients
                    term = J_kron_Jconj[k, l]
                    if np.isclose(coeff, 1):
                        expr.append(f"+{term}")
                    elif np.isclose(coeff, -1):
                        expr.append(f"-{term}")
                    elif np.isclose(coeff, 1j):
                        expr.append(f"+I*{term}")
                    elif np.isclose(coeff, -1j):
                        expr.append(f"-I*{term}")
        # Ensure expr is not empty
        M[i, j] = "".join(expr).lstrip('+').replace('+-', '-') if expr else "Empty"

# Debug: Print raw M to see what's being computed
print("Raw Mueller Matrix M (before formatting):")
for i in range(4):
    for j in range(4):
        print(f"M_{i+1}{j+1} = {M[i, j]}")
    print()

# Format the Mueller matrix elements
def format_element(expr, i, j):
    if not expr or expr == "Empty":
        return "Error: Empty expression"
    
    terms = [t.strip() for t in expr.replace('-', '+-').split('+') if t.strip()]
    formatted_terms = []
    
    for term in terms:
        if 'I*' in term:
            formatted_terms.append(term.replace('I*', ''))
        elif '*' in term:
            a, b = term.split('*')
            if a == b.replace('*', ''):
                formatted_terms.append(f"|{a}|^2")
            else:
                formatted_terms.append(term)
        else:
            formatted_terms.append(term)
    
    result = ' + '.join(formatted_terms).replace(' + -', ' - ')
    
    # Apply Re or Im based on position
    if (i, j) in [(0, 0), (1, 1), (3, 3)]:
        return result
    elif (i, j) in [(0, 3), (1, 3), (2, 3), (3, 0), (3, 1), (3, 2)]:
        return f"Im({result})"
    else:
        return f"Re({result})"

# Display the formatted Mueller matrix
print("Formatted Mueller Matrix M:")
for i in range(4):
    for j in range(4):
        formatted = format_element(M[i, j], i, j)
        print(f"M_{i+1}{j+1} = {formatted}")
    print()