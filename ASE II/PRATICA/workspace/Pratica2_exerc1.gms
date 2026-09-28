$Title Despacho Economico - Parte 1

Sets
    i geradores / G1, G2, G3 /;

Parameters
    a(i) coeficiente quadratico
    / G1 0.004, G2 0.006, G3 0.008 /
    
    b(i) coeficiente linear
    / G1 5.2, G2 5.0, G3 5.0 /
    
    c(i) custo fixo
    / G1 500, G2 400, G3 200 /
    
    Pmin(i) limite minimo (MW)
    / G1 200, G2 150, G3 100 /
    
    Pmax(i) limite maximo (MW)
    / G1 450, G2 350, G3 235 /
    
    Pd demanda do sistema (MW)
    / 500 /;

Variables
    Pg(i) potencia gerada (MW)
    CT    custo total do sistema (dolares por hora);

* Define que a geracao nao pode ser negativa
Positive Variable Pg;

Equations
    Funcao_Custo   equacao da funcao objetivo
    Balanco_Pot    restricao de atendimento da demanda;

* Funcao objetivo: Somatorio dos custos individuais
Funcao_Custo..
    CT =e= sum(i, c(i) + b(i)*Pg(i) + a(i)*power(Pg(i),2));

* Balanco de Potencia: A soma da geracao deve ser igual a demanda
Balanco_Pot..
    sum(i, Pg(i)) =e= Pd;

* Aplicacao dos limites operacionais maximos e minimos
Pg.lo(i) = Pmin(i);
Pg.up(i) = Pmax(i);

* Definicao e resolucao do modelo
Model Despacho /all/;

* NLP (Non-Linear Programming) e utilizado devido ao termo quadratico na funcao de custo
Solve Despacho using nlp minimizing CT;

* Exibe os resultados no arquivo .lst
Display Pg.l, CT.l;
