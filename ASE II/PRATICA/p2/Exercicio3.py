import numpy as np

# Dados dos Geradores
a_cost = np.array([0.004, 0.006, 0.008]) # Coeficientes quadráticos originais
b_cost = np.array([5.2, 5.0, 5.0])       # Coeficientes lineares
c_cost = np.array([500.0, 400.0, 200.0]) # Custos fixos

# O Slide 31 define dF/dP = a*P + b. 
# Logo, o 'a' da matriz do slide equivale a 2*a_cost.
a_slide = 2 * a_cost 
b_slide = b_cost

# Matriz B do caso base (em p.u. - Base 100 MVA)
Sbase = 100.0
B_pu = np.array([[ 0.00693144, -0.00004987, -0.00011201],
                 [-0.00004987,  0.00951274,  0.00081946],
                 [-0.00011201,  0.00081946,  0.00344373]])
B0_pu = np.array([0.00060614, 0.00098124, -0.00010308])
B00_pu = 0.00009021

# Conversão da Matriz B para operar diretamente em MW (unidade do despacho)
B = B_pu / Sbase
B0 = B0_pu
B00 = B00_pu * Sbase

# PASSO 1: Especificar o nível de carga do sistema
Pd = 500.0 
tolerancia = 1e-4

# PASSO 2: Primeira aproximação de lambda (supondo perdas zero)
# Equação sem perdas: sum((lambda - b) / a_slide) = Pd
lambda_k = (Pd + np.sum(b_slide / a_slide)) / np.sum(1 / a_slide)
lambda_k = 6.938462

# Variável auxiliar para o método da secante (Atualização de lambda do Slide 32)
lambda_k_minus_1 = lambda_k - 0.1 
P_k_minus_1 = np.zeros(3)

print("Iniciando processo iterativo...\n")

for k in range(1, 50): # Limite de 50 iterações para segurança
    # PASSO 3: Resolver o sistema de equações para obter P_gi
    # Montagem da Matriz KxK e Vetor Solução baseada no Slide 32
    Matriz_A = np.zeros((3, 3))
    Vetor_C = np.zeros(3)
    
    for i in range(3):
        for j in range(3):
            if i == j:
                # Diagonal: (a_i / lambda) + 2*B_ii
                Matriz_A[i,j] = (a_slide[i] / lambda_k) + 2 * B[i,j] 
            else:
                # Fora da diagonal: 2*B_ij
                Matriz_A[i,j] = 2 * B[i,j] 
        
        # Vetor do lado direito: (1 - B_i0) - (b_i / lambda)
        Vetor_C[i] = (1 - B0[i]) - (b_slide[i] / lambda_k) 

    # Resolve o sistema linear A * P = C
    P_k = np.linalg.solve(Matriz_A, Vetor_C)
    
    # PASSO 4: Calcular as perdas de transmissão
    # PL = P^T * B * P + B0^T * P + B00
    PL_k = P_k.T @ B @ P_k + B0.T @ P_k + B00 
    
    # PASSO 5: Verificar o balanço de potência
    erro = Pd + PL_k - np.sum(P_k) 
    
    print(f"Iteração {k}: Lambda = {lambda_k:.5f}, Erro = {erro:.5f} MW")
    
    # Se o balanço ocorrer dentro da tolerância, finaliza o processo
    if abs(erro) < tolerancia:
        print("\n=== CONVERGÊNCIA ALCANÇADA ===")
        break
        
    # Atualizar o valor de lambda (Fórmula exata do Slide 32 - Método da Secante)
    if k > 1:
        numerador = lambda_k - lambda_k_minus_1
        denominador = np.sum(P_k) - np.sum(P_k_minus_1)
        delta_lambda = (numerador / denominador) * erro 
    else:
        # Pulo inicial para viabilizar a primeira diferença do secante
        delta_lambda = 0.05 * erro 
        
    lambda_k_minus_1 = lambda_k
    P_k_minus_1 = P_k.copy()
    
    # PASSO 6: Retorna para a etapa 3
    lambda_k = lambda_k + delta_lambda 

# Resultados Finais
print(f"Geração P1 = {P_k[0]:.2f} MW")
print(f"Geração P2 = {P_k[1]:.2f} MW")
print(f"Geração P3 = {P_k[2]:.2f} MW")
print(f"Perdas PL  = {PL_k:.4f} MW")
print(f"Geração Total = {np.sum(P_k):.2f} MW (Demanda + Perdas)")
print(lambda_k)

# Custo Total
CT = np.sum(c_cost + b_cost*P_k + a_cost*(P_k**2))
print(f"Custo Total = ${CT:.2f}/h")

# === CÁLCULO DOS FATORES DE PENALIZAÇÃO E CUSTOS INCREMENTAIS ===
print("\n=== ANÁLISE DE CUSTOS E PENALIZAÇÃO ===")
for i in range(3):
    # Custo Incremental: dF/dP = 2*a*P + b (Lembrando que a_slide já é 2*a)
    custo_inc = a_slide[i] * P_k[i] + b_slide[i]
    
    # Fator de Penalização: L = lambda / Custo Incremental
    fator_penalizacao = lambda_k / custo_inc
    
    print(f"Gerador {i+1}:")
    print(f"  Custo Incremental (IC) = ${custo_inc:.4f}/MWh")
    print(f"  Fator de Penalização (L) = {fator_penalizacao:.4f}")