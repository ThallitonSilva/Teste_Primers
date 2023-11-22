import streamlit as st
import streamlit_ext as ste
from test_funcs import *

st.set_page_config(page_title='Teste dos Primers', layout='wide')

hide_menu_style = """
        <style>
        #MainMenu {visibility: hidden; }
        footer {visibility: hidden;}
        </style>
        """

st.markdown(hide_menu_style, unsafe_allow_html=True)

title_col1, title_col2, title_col3, title_col4 = st.columns((1, 2, 2, 1))
col1, col2, col3, col4 = st.columns((1, 2, 2, 1))

title_col2.title('Teste dos primers')

primers = col2.file_uploader('Insira um arquivo excel com os primers para serem testados',
                             type='xlsx',
                             accept_multiple_files=False)

genome = col3.file_uploader('Insira o genoma (formato FASTA) no qual os primers serão testados',
                            type='fasta',
                            accept_multiple_files=False)

exp = col2.expander("A tabela deve, obrigatoriamente, ter as colunas Name, Primer_ID, Forward e Reverse. Exemplo:")

exemplo = [['Name', 'Primer_ID', 'Forward', 'Reverse'],
           ['Primer_1', 'AT1G60590', 'CGCAGCTACAACCCACACTA', 'TGGTGGGTTCTCAACATCAA'],
           ['Primer_2', 'AT4G37990', 'TGCCATCCTTGTCTTGTACG', 'CCGCTCCAAACACTCTTCTC'],
           ['Primer_3', 'AT2G33380', 'GAGGCACAAATCCCTCAAGC', 'CCGAAAAGTACCAATATGTAACG']]
ex = pd.DataFrame(exemplo[1:], columns=exemplo[0])
exp.write(ex)

if primers and genome:

    table_primers = pd.read_excel(primers)

    col2.write(table_primers)

    list_fasta = fasta_genes(genome)

    col3.write(f'O arquivo fasta inserido possui {len(list_fasta)} sequencias')

    if ('Primer_ID' and 'Forward' and 'Reverse') not in table_primers.columns:
        st.write('O nome das colunas estão errados. Por favor, arrume!')
        st.write('Veja o exemplo fornecido.')
        col2.write(table_primers)

    else:

        st.write('Iniciando o Teste dos primers...')

        ok, nao_ok, amp, amp_all = verifica_primers(table_primers, list_fasta)

        passou = faz_tabela(ok)
        n_passou = faz_tabela(nao_ok)

        final = pd.concat([passou, n_passou]).sort_values(by=['Contagem', 'Gene_ID'])

        st.write('Primers Testados...')

        st.write(final)

        excel = make_excel(final)

        ste.download_button(label="Download Resultados",
                            data=excel,
                            file_name="Teste_dos_Primers.xlsx",
                            mime="application/vnd.ms-excel")

        arq = desenha_gel(amp, amp_all)

        ste.download_button(label="Download Imagens",
                            data=arq,
                            file_name="Imagens.zip")

