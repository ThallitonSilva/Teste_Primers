import pandas as pd
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import Anneal
from pydna.readers import read
from pydna.gel import gel, interpolator
from pydna.ladders import HI_LO_DNA_MARKER
import streamlit as st
import io
import zipfile


def fasta_genes(fasta):
    stringio = io.StringIO(fasta.getvalue().decode("utf-8"))

    seqs = SeqIO.parse(stringio, 'fasta')

    return [seq.upper() for seq in seqs]


def pcr_in_silico(fasta, table, limit):
    temp = []

    for j in table.iloc:

        for gene in fasta:

            if j['ID'] == gene.id:
                sequence = Dseqrecord(gene.seq)

                fwd = read(f'>{j.ID}\n{j.New_Forward}', ds=False)
                rev = read(f'>{j.ID}\n{j.New_Reverse}', ds=False)

                pcr = Anneal((fwd, rev), sequence, limit=limit)

                if len(pcr.products) == 1:
                    temp.append(j)

    return pd.DataFrame(temp)


def verifica_primers(dataframe_primers, fasta):
    passou_teste = []
    nao_passou = []

    temp_amplicons_un = []
    temp_amplicons_all = []

    progress_text = "Testando Primers..."
    my_bar = st.progress(0, text=progress_text)

    tamanho_fragmento = (100 / dataframe_primers.shape[0]) / 100
    percent = 0

    for primer in dataframe_primers.iloc:

        temp_amplicons = []

        fwd = read(f'>{primer.Primer_ID}_F\n{primer.Forward}', ds=False)
        rev = read(f'>{primer.Primer_ID}_R\n{primer.Reverse}', ds=False)

        for gene in fasta:

            sequence = Dseqrecord(gene.seq)

            pcr_silico = Anneal((fwd, rev), sequence)

            if len(pcr_silico.products) >= 1:
                temp_amplicons.append([gene.id, pcr_silico])
                temp_amplicons_un.append([primer.Primer_ID, pcr_silico.products])
                temp_amplicons_all.append(pcr_silico.products)

        count = len(temp_amplicons)

        if count == 0:
            passou_teste.append([primer.Primer_ID, count, f'O primer {primer.Primer_ID} não se anelou em nenhum gene'])

        if count == 1:

            for ampli in temp_amplicons:
                gene = ampli[0]
                texto = str(ampli[1])
                #p1, p2, p3 = texto.split('\n')
                #fw = p2.split()[-1]
                #rv = p3.split()[-1]

            passou_teste.append([primer.Primer_ID, count, f'O primer {primer.Primer_ID} se anelou somente no gene {gene}'])

        if count > 1:

            opcs = []

            for ampli in temp_amplicons:
                gene = ampli[0]
                texto = str(ampli[1])
                #p1, p2, p3 = texto.split('\n')
                #fw = p2.split()[-1]
                #rv = p3.split()[-1]
                opcs.append(f'O primer {primer.Primer_ID} se anelou no gene {gene}')

            nao_passou.append([primer.Primer_ID, count, ' ;'.join(opcs)])

        percent += tamanho_fragmento
        if percent >= 1:
            percent = 0.99
        
        my_bar.progress(percent, text=f'Primer {primer.Primer_ID} testado!')

    return passou_teste, nao_passou, temp_amplicons_un, temp_amplicons_all


# Crio um DataFrame a partir da lista de primers que passaram, ou não, no teste
def faz_tabela(lista_teste):
    return pd.DataFrame(lista_teste,
                        columns=['Gene_ID', 'Contagem', 'Descrição'])


# Desenha uma simulação dos geis de agarose, com base nos primers e nos genes
# Caso um primer se anele em mais de um gene, o gel terá várias bandas, uma para cada anelamento
def desenha_gel(temp_amplicons, temp_amplicons_all):
    buff = io.BytesIO()

    with zipfile.ZipFile(buff, 'w') as zp:

        buf = io.BytesIO()

        try:
            for temp in temp_amplicons:
                t1, *t2 = temp
                gel_im = gel([HI_LO_DNA_MARKER, *t2], interpolator=interpolator(HI_LO_DNA_MARKER))
                gel_im.save(fp=buf, format='png')
                buf.seek(0)
                zp.writestr(f'Gel_{t1}.png', buf.getvalue())

        except:
            msg = f'Não foi possível desenhar o gel'
            st.write(msg)

        try:
            gel_im = gel([HI_LO_DNA_MARKER, *temp_amplicons_all], interpolator=interpolator(HI_LO_DNA_MARKER))
            gel_im.save(fp=buf, format='png')
            buf.seek(0)
            zp.writestr('Gel_All.png', buf.getvalue())

        except:
            st.write('Não foi possível desenhar o gel de todos os genes')

    return buff

def make_excel(table):

    buffer = io.BytesIO()

    with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:

        table.to_excel(writer, sheet_name=f'Teste_Primers', index=False)

    #writer.save()

    return buffer
