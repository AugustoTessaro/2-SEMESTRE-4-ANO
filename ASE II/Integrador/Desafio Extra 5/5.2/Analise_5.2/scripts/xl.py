import json,os
from openpyxl import Workbook
from openpyxl.styles import Font,PatternFill,Alignment,Border,Side
from openpyxl.utils import get_column_letter as L
from openpyxl.formatting.rule import FormulaRule
R=json.load(open(os.path.expanduser('~/work/crit.json')))
OUT=os.path.expanduser('~/mnt/5.2/Analise_5.2/Item_5.2_Criterios_Estabilidade.xlsx')
wb=Workbook()
F=Font(name='Arial',size=10); FB=Font(name='Arial',size=10,bold=True); FBLUE=Font(name='Arial',size=10,color='0000FF')
HDR=PatternFill('solid',fgColor='DDE6F0'); thin=Side(style='thin',color='BFBFBF'); BD=Border(top=thin,bottom=thin,left=thin,right=thin)
RED=PatternFill('solid',fgColor='F8CBAD'); GRN=PatternFill('solid',fgColor='E2EFDA')
EVN={'E1_ABRE_LT9-14':'E1 – abertura LT 9–14','E2_FT_LT2-4_150ms':'E2 – curto 1F LT 2–4 (150 ms)'}
TE={'E1_ABRE_LT9-14':1.0,'E2_FT_LT2-4_150ms':1.15}
def hdr(ws,row,cols):
    cols=cols if isinstance(cols,list) else [ws.cell(row,j).value for j in range(1,cols+1)]
    for j,c in enumerate(cols,1):
        x=ws.cell(row,j,c); x.font=FB; x.fill=HDR; x.border=BD; x.alignment=Alignment(wrap_text=True,vertical='center',horizontal='center')
    ws.row_dimensions[row].height=42
def style(ws,r0,r1,ncol):
    for r in range(r0,r1+1):
        for c in range(1,ncol+1):
            x=ws.cell(r,c); x.border=BD
            if x.font!=FBLUE: x.font=F if not (x.font and x.font.color and x.font.color.rgb=='FF0000FF') else x.font
# ---------- Premissas
P=wb.active; P.title='Premissas'
rows=[('Item','Valor','Referência / observação'),
('Tensão mínima na 1ª oscilação pós-distúrbio (pu)',0.60,'Sub. 2.3, 2.4.3.7 (a) – nenhuma barra é de 500 kV'),
('Tensão mínima nas demais oscilações (pu)',0.80,'Sub. 2.3, 2.4.3.7 (b)'),
('Máxima queda de tensão entre início e fim (pu)',0.10,'Sub. 2.3, 2.4.3.7 (c): Vfinal ≥ Vinicial − 10 % Vnop'),
('Amplitude máxima pico a pico após 10 s da eliminação (pu)',0.02,'Sub. 2.3, 2.4.3.7 (d)'),
('Faixa de regime final – mínimo (pu)',0.95,'Sub. 2.3, Tabela 1 (DGLT 2, condição normal); emergência admite 0,90'),
('Faixa de regime final – máximo (pu)',1.05,'Sub. 2.3, Tabela 1 (DGLT 2)'),
('1º estágio do ERAC (Hz)',59.50,'Arquivo .stb fornecido: 0,9917 × 60 Hz'),
('Tempo de eliminação E1 (s)',1.00,'Abertura sem falta em t = 1 s'),
('Tempo de eliminação E2 (s)',1.15,'Curto em 1 s + 150 ms (Sub. 2.3, Tabela 2, 138 kV)'),
('Fim da simulação (s)',20.0,'Janela do critério (d): de t_elim + 10 s até 20 s'),
('Definição da 1ª oscilação','—','Da eliminação até o 2º extremo do ângulo do gerador mais perturbado (E1 antes: 1,00–1,75 s; E1 depois: 1,00–1,91 s; E2 antes: 1,15–1,57 s; E2 depois: 1,15–1,42 s)'),
('Barras analisadas','—','Todas as físicas (1–6, 8–14); barra 7 é fictícia (ponto central do trafo de 3 enrolamentos)'),
('Limites das LTs','—','DLIN do PWF: 50 MVA (66 kV) e 100 MVA (138 kV); condição normal = emergência (sem valor distinto)'),
('Fonte dos valores numéricos','—','Arquivos .plt do ANATEM v13 (pasta 5.2.1), passo de plotagem 5 ms; extraídos por script'),
]
for i,r in enumerate(rows,1):
    for j,v in enumerate(r,1): P.cell(i,j,v)
hdr(P,1,3)
for i in range(2,len(rows)+1):
    for j in range(1,4): P.cell(i,j).font=F; P.cell(i,j).border=BD; P.cell(i,j).alignment=Alignment(wrap_text=True,vertical='top')
    if isinstance(rows[i-1][1],float): P.cell(i,2).font=FBLUE
P.column_dimensions['A'].width=48; P.column_dimensions['B'].width=10; P.column_dimensions['C'].width=80
# names
Vmin1='Premissas!$B$2'; Vmin2='Premissas!$B$3'; dVlim='Premissas!$B$4'; pplim='Premissas!$B$5'; Vlo='Premissas!$B$6'; Vhi='Premissas!$B$7'
# ---------- Tensoes
T=wb.create_sheet('Tensoes')
cols=['Evento','Caso','Barra','Tipo','V inicial (pu)','V mín. durante a falta (pu)','V mín. 1ª oscilação (pu)','V mín. demais oscilações (pu)','V final (pu)','V máx. pós-evento (pu)','Pico a pico após t_elim+10 s (pu)',
      '(a) ≥ 0,60','(b) ≥ 0,80','(c) Vfin ≥ Vini − 0,10','(d) pp ≤ 0,02','Vfinal em 0,95–1,05']
hdr(T,1,cols); r=2
GEN={1,2,3,6,8,12,14}
for ev in EVN:
  for c in ['ANTES','DEPOIS']:
    for b,v in R[f'{ev}|{c}']['V'].items():
      b=int(b)
      vals=[EVN[ev],c.capitalize(),b,'Geração' if b in GEN else 'Carga',v['Vini'],v['Vfalta_min'] if v['Vfalta_min'] else '—',v['Vmin_1a'],v['Vmin_dem'],v['Vfin'],v['Vmax_pos'],v['pp10']]
      for j,x in enumerate(vals,1): T.cell(r,j,x)
      T.cell(r,12,f'=IF(G{r}>={Vmin1},"OK","VIOLA")'); T.cell(r,13,f'=IF(H{r}>={Vmin2},"OK","VIOLA")')
      T.cell(r,14,f'=IF(I{r}>=E{r}-{dVlim},"OK","VIOLA")'); T.cell(r,15,f'=IF(K{r}<={pplim},"OK","VIOLA")')
      T.cell(r,16,f'=IF(AND(ROUND(I{r},2)>={Vlo},ROUND(I{r},2)<={Vhi}),"OK","FORA")')
      for j in range(5,12): T.cell(r,j).number_format='0.0000' if j!=11 else '0.00000'
      r+=1
TN=r-1
for rr in range(2,TN+1):
    for j in range(1,17): x=T.cell(rr,j); x.font=F; x.border=BD; x.alignment=Alignment(horizontal='center')
T.conditional_formatting.add(f'L2:P{TN}',FormulaRule(formula=['OR(L2="VIOLA",L2="FORA")'],fill=RED))
T.conditional_formatting.add(f'L2:P{TN}',FormulaRule(formula=['L2="OK"'],fill=GRN))
for j,w in enumerate([26,9,7,9]+[11]*12,1): T.column_dimensions[L(j)].width=w
T.cell(TN+2,1,'Faixa 0,95–1,05 verificada com 2 casas decimais (precisão da Tabela 1). Única ocorrência acima de 1,050: barra 14 no E1 antes, 1,0506 pu (erro de regime do regulador da própria barra, setpoint 1,05).').font=Font(name='Arial',size=9,italic=True)
T.freeze_panes='E2'
# ---------- Angulos
A=wb.create_sheet('Angulos')
cols=['Evento','Caso','Gerador','δ−δ1 inicial (°)','δ−δ1 mín. pós (°)','δ−δ1 máx. pós (°)','δ−δ1 final (°)','Excursão máx. a partir do inicial (°)','Deslocamento final (°)','Pico a pico últimos 5 s (°)','Sincronismo']
hdr(A,1,cols); r=2
for ev in EVN:
  for c in ['ANTES','DEPOIS']:
    q=R[f'{ev}|{c}']
    for g in ['2','3','6','8','12','14']:
      for j,x in enumerate([EVN[ev],c.capitalize(),int(g),q['ang_ini'][g],q['ang_min'][g],q['ang_max'][g],q['ang_fin'][g]],1): A.cell(r,j,x)
      A.cell(r,8,f'=MAX(F{r}-D{r},D{r}-E{r})'); A.cell(r,9,f'=G{r}-D{r}'); A.cell(r,10,q['ang_pp_ult5s'][g])
      A.cell(r,11,f'=IF(J{r}<0.5,"Mantido","VERIFICAR")')
      for j in range(4,11): A.cell(r,j).number_format='0.00' if j!=10 else '0.000'
      r+=1
AN=r-1
for rr in range(2,AN+1):
    for j in range(1,12): x=A.cell(rr,j); x.font=F; x.border=BD; x.alignment=Alignment(horizontal='center')
A.cell(AN+2,1,'Sincronismo: oscilações amortecidas, com pico a pico < 0,5° nos últimos 5 s (critério adotado para caracterizar retorno a regime). Referência comum: máquina da barra 1.').font=Font(name='Arial',size=9,italic=True)
for j,w in enumerate([26,9,9]+[13]*8,1): A.column_dimensions[L(j)].width=w
# ---------- Frequencia
Q=wb.create_sheet('Frequencia')
cols=['Evento','Caso','Gerador','f mín. pós (Hz)','f máx. pós (Hz)','f final (Hz)','Margem ao 1º estágio do ERAC (Hz)']
hdr(Q,1,cols); r=2
for ev in EVN:
  for c in ['ANTES','DEPOIS']:
    q=R[f'{ev}|{c}']
    for g in ['1','2','3','6','8','12','14']:
      for j,x in enumerate([EVN[ev],c.capitalize(),int(g),q['f_min'][g],q['f_max'][g],q['f_fin'][g]],1): Q.cell(r,j,x)
      Q.cell(r,7,f'=D{r}-Premissas!$B$8')
      for j in range(4,8): Q.cell(r,j).number_format='0.000'
      r+=1
QN=r-1
for rr in range(2,QN+1):
    for j in range(1,8): x=Q.cell(rr,j); x.font=F; x.border=BD; x.alignment=Alignment(horizontal='center')
for j,w in enumerate([26,9,9,13,13,13,16],1): Q.column_dimensions[L(j)].width=w
# ---------- Fluxos
X=wb.create_sheet('Fluxos')
cols=['Evento','Caso','LT','S inicial (MVA)','S máx. pós-evento (MVA)','S final (MVA)','Limite (MVA)','Carregamento inicial (%)','Carregamento final (%)','Situação final']
hdr(X,1,cols); r=2
for ev in EVN:
  for c in ['ANTES','DEPOIS']:
    for k,f in R[f'{ev}|{c}']['fluxo'].items():
      for j,x in enumerate([EVN[ev],c.capitalize(),k,f['S0'],f['Smax_pos'],f['Sfin']],1): X.cell(r,j,x)
      X.cell(r,7,f['lim']).font=FBLUE
      X.cell(r,8,f'=D{r}/G{r}'); X.cell(r,9,f'=F{r}/G{r}')
      X.cell(r,10,f'=IF(F{r}=0,"Aberta",IF(I{r}>1,"SOBRECARGA","OK"))')
      for j in (4,5,6): X.cell(r,j).number_format='0.0'
      for j in (8,9): X.cell(r,j).number_format='0.0%'
      r+=1
XN=r-1
for rr in range(2,XN+1):
    for j in range(1,11):
        x=X.cell(rr,j); x.border=BD; x.alignment=Alignment(horizontal='center')
        if j!=7: x.font=F
X.conditional_formatting.add(f'J2:J{XN}',FormulaRule(formula=['J2="SOBRECARGA"'],fill=RED))
for j,w in enumerate([26,9,8,11,13,11,10,13,13,13],1): X.column_dimensions[L(j)].width=w
# ---------- Resumo
S=wb.create_sheet('Resumo',0)
S.cell(1,1,'Item 5.2 – Estudo de estabilidade eletromecânica: verificação dos critérios (Sub. 2.3, 2.4.3)').font=Font(name='Arial',size=12,bold=True)
cols=['Evento','Caso','Sincronismo','Maior excursão angular (°)','f mínima (Hz)','V mín. 1ª osc. (pu)','(a)','V mín. demais (pu)','(b)','Maior queda Vini−Vfin (pu)','(c)','Maior pico a pico (pu)','(d)','Barras fora de 0,95–1,05','Maior carregamento final (%)','LTs em sobrecarga']
hdr(S,3,cols); r=4
for ev in EVN:
  for c in ['ANTES','DEPOIS']:
    ce=c.capitalize(); e=EVN[ev]
    S.cell(r,1,e); S.cell(r,2,ce)
    cond=f'Angulos!$A$2:$A${AN},$A{r},Angulos!$B$2:$B${AN},$B{r}'
    S.cell(r,3,f'=IF(COUNTIFS({cond},Angulos!$K$2:$K${AN},"Mantido")=6,"Mantido","VERIFICAR")')
    S.cell(r,4,f'=_xlfn.MAXIFS(Angulos!$H$2:$H${AN},{cond})')
    S.cell(r,5,f'=_xlfn.MINIFS(Frequencia!$D$2:$D${QN},Frequencia!$A$2:$A${QN},$A{r},Frequencia!$B$2:$B${QN},$B{r})')
    tc=f'Tensoes!$A$2:$A${TN},$A{r},Tensoes!$B$2:$B${TN},$B{r}'
    S.cell(r,6,f'=_xlfn.MINIFS(Tensoes!$G$2:$G${TN},{tc})'); S.cell(r,7,f'=IF(F{r}>={Vmin1},"OK","VIOLA")')
    S.cell(r,8,f'=_xlfn.MINIFS(Tensoes!$H$2:$H${TN},{tc})'); S.cell(r,9,f'=IF(H{r}>={Vmin2},"OK","VIOLA")')
    S.cell(r,10,f'=SUMPRODUCT(MAX((Tensoes!$A$2:$A${TN}=$A{r})*(Tensoes!$B$2:$B${TN}=$B{r})*(Tensoes!$E$2:$E${TN}-Tensoes!$I$2:$I${TN})))'); S.cell(r,11,f'=IF(J{r}<={dVlim},"OK","VIOLA")')
    S.cell(r,12,f'=_xlfn.MAXIFS(Tensoes!$K$2:$K${TN},{tc})'); S.cell(r,13,f'=IF(L{r}<={pplim},"OK","VIOLA")')
    S.cell(r,14,f'=COUNTIFS({tc},Tensoes!$P$2:$P${TN},"FORA")')
    xc=f'Fluxos!$A$2:$A${XN},$A{r},Fluxos!$B$2:$B${XN},$B{r}'
    S.cell(r,15,f'=_xlfn.MAXIFS(Fluxos!$I$2:$I${XN},{xc})')
    S.cell(r,16,f'=COUNTIFS({xc},Fluxos!$J$2:$J${XN},"SOBRECARGA")')
    for j,fmt in {4:'0.0',5:'0.000',6:'0.000',8:'0.000',10:'0.000',12:'0.0000',15:'0.0%'}.items(): S.cell(r,j).number_format=fmt
    r+=1
for rr in range(4,r):
    for j in range(1,17): x=S.cell(rr,j); x.font=F; x.border=BD; x.alignment=Alignment(horizontal='center',vertical='center')
S.conditional_formatting.add(f'C4:P{r-1}',FormulaRule(formula=['OR(C4="VIOLA",C4="VERIFICAR")'],fill=RED))
S.conditional_formatting.add(f'O4:O{r-1}',FormulaRule(formula=['O4>1'],fill=RED))
S.conditional_formatting.add(f'P4:P{r-1}',FormulaRule(formula=['P4>0'],fill=RED))
for j,w in enumerate([28,9,12,12,10,11,6,11,6,13,6,11,6,12,13,11],1): S.column_dimensions[L(j)].width=w
notas=['Notas:',
'• E1 depois: sincronismo e tensões atendem, mas a LT 13–14 termina com 83 MVA (166 %) e a LT 6–13 com 61,6 MVA (123 %); viola 2.4.3.4/2.4.3.5. Aumento de perdas (~8 MW) leva a frequência a 59,71 Hz (G14) e 59,93 Hz em 20 s, sem atingir o ERAC.',
'• E2 antes: a LT 4–5 termina com 106,5 MVA (107 %); no depois, 92 %. A conexão reduz a excursão angular (G12: 9,1° → 5,0°) e eleva as tensões na falta e na recuperação.',
'• Tensão durante a falta (E2) não entra nos critérios (a)–(d), que são pós-distúrbio; a menor é 0,73 pu (antes) e 0,75 pu (depois) na barra 4.']
for i,n in enumerate(notas): S.cell(r+1+i,1,n).font=Font(name='Arial',size=9,bold=(i==0))
S.freeze_panes='C4'
wb.save(OUT); print(OUT)
