import os
import re
import requests
import pandas as pd
from Bio import Entrez, SeqIO
from typing import Dict

class NCBIGenomeMiner:
    """
    Classe para mineração de dados genômicos utilizando a biblioteca Bio.Entrez.
    
    Esta classe encapsula a comunicação oficial com a API Entrez do NCBI para:
    1. Mapear cromossomos e metadados de montagens genômicas de referência.
    2. Localizar o locus e coordenadas genômicas de genes específicos.
    3. Fazer download de sequências cromossômicas em formato FASTA de maneira segura.
    """

    def __init__(self, taxon: str):
        """
        Inicializa o minerador genômico para uma espécie específica.
        
        Parameters:
            taxon (str): Nome científico da espécie (ex: 'Gallus gallus').
        """
        self.taxon = taxon
        
        # [BOAS PRÁTICAS NCBI] O NCBI exige e-mail e identificação da ferramenta
        # para evitar bloqueio do endereço IP por acessos automatizados descontrolados.
        Entrez.email = "albermangomes3@gmail.com"
        Entrez.tool = "ProjetoBioinfo2026"

    def get_chromosome_metadata(self) -> pd.DataFrame:
        """
        [PASSO 1 DO PIPELINE]
        Obtém os metadados dos cromossomos da espécie a partir do genoma de referência.
        
        Returns:
            pd.DataFrame: Tabela (Tabela 1) contendo Espécie, Nome do Cromossomo,
                          Tamanho em pares de bases (bp) e ID do RefSeq (FASTA).
        """
        print(f"--- [PASSO 1] Consultando genoma de {self.taxon} via Entrez ---")
        chromosomes_data = []

        try:
            # -------------------------------------------------------------
            # Etapa 1.1: Busca o ID da montagem de referência (Assembly)
            # -------------------------------------------------------------
            print("1. Buscando o genoma de referência (Assembly)...")
            search_handle = Entrez.esearch(
                db="assembly", 
                term=f'"{self.taxon}"[Organism] AND "reference genome"[filter]', 
                retmax=1
            )
            assembly_record = Entrez.read(search_handle)
            search_handle.close()

            # Valida se o NCBI retornou alguma montagem válida para a espécie
            if not assembly_record.get("IdList"):
                raise ValueError(f"Nenhum genoma de referência encontrado para '{self.taxon}'.")

            # -------------------------------------------------------------
            # Etapa 1.2: Obtém o código de acesso da montagem (ex: GCF_016699485.2)
            # -------------------------------------------------------------
            sum_handle = Entrez.esummary(db="assembly", id=assembly_record["IdList"][0])
            assembly_summary = Entrez.read(sum_handle)
            sum_handle.close()
            
            assembly_acc = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            print(f"Montagem encontrada com sucesso: {assembly_acc}")

            # -------------------------------------------------------------
            # Etapa 1.3: Mapeia as sequências genômicas associadas à montagem
            # -------------------------------------------------------------
            print("2. Mapeando os cromossomos (Nuccore)...")
            nuc_search = Entrez.esearch(
                db="nuccore", 
                term=f'{assembly_acc}[Assembly] AND biomol_genomic[PROP]', 
                retmax=500
            )
            nuc_record = Entrez.read(nuc_search)
            nuc_search.close()
            
            seq_ids = nuc_record.get("IdList", [])
            
            # -------------------------------------------------------------
            # Etapa 1.4: Extrai os detalhes de cada sequência (tamanho e título)
            # -------------------------------------------------------------
            nuc_sum = Entrez.esummary(db="nuccore", id=",".join(seq_ids))
            seq_summaries = Entrez.read(nuc_sum)
            nuc_sum.close()

            for seq in seq_summaries:
                acc = seq.get("AccessionVersion", "")
                length = seq.get("Length", 0)
                title = seq.get("Title", "")
                
                # Filtra apenas sequências válidas maiores que 1 bp
                if length and int(length) > 1:
                    # Aplica expressão regular no título para identificar se é cromossomo
                    if "chromosome" in title.lower():
                        match_chr = re.search(r'(chromosome\s+\w+)', title, re.IGNORECASE)
                        chr_name = match_chr.group(1) if match_chr else "Cromossomo"
                    else:
                        chr_name = "Scaffold"
                        
                    chromosomes_data.append({
                        "Espécie": self.taxon,
                        "Cromossomo": chr_name,
                        "Tamanho (bp)": int(length),
                        "ID do Cromossomo (FASTA)": acc
                    })

        except Exception as e:
            print(f"\n[ERRO NA API] Detalhes: {e}")
            print("Ocorreu uma falha de comunicação com o servidor Entrez.")

        # Converte a lista em DataFrame do Pandas
        df = pd.DataFrame(chromosomes_data)
        if df.empty:
            raise RuntimeError(f"Falha definitiva ao extrair cromossomos para {self.taxon}.")
            
        # Remove duplicatas de ID de acesso
        df = df.drop_duplicates(subset=["ID do Cromossomo (FASTA)"])

        # Filtro estrito: Mantém apenas cromossomos completos (elimina contigs e scaffolds)
        df = df[df["Cromossomo"].str.contains("chromosome", case=False, na=False)]
        
        # Ordena a tabela do maior para o menor cromossomo
        return df.sort_values(by="Tamanho (bp)", ascending=False).reset_index(drop=True)

    def locate_gene(self, gene_symbol: str) -> Dict[str, str]:
        """
        [PASSO 2 DO PIPELINE]
        Localiza o gene no banco 'gene' do NCBI e recupera seu locus genômico
        (ID do cromossomo, número do cromossomo e posições de início e fim).
        
        Parameters:
            gene_symbol (str): Símbolo do gene (ex: 'GAPDH').
            
        Returns:
            Dict[str, str]: Dicionário contendo os metadados e coordenadas do gene.
        """
        print(f"\n--- [PASSO 2] Localizando coordenadas do gene {gene_symbol} ---")
        
        # Busca o ID do gene no banco de dados 'gene' do NCBI
        search_handle = Entrez.esearch(
            db="gene", 
            term=f'"{self.taxon}"[Organism] AND {gene_symbol}[Gene Name]', 
            retmax=1
        )
        record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = record.get("IdList", [])
        if not id_list:
            raise ValueError(f"Gene '{gene_symbol}' não encontrado para {self.taxon}.")
            
        # Consulta o sumário detalhado do gene encontrado
        sum_handle = Entrez.esummary(db="gene", id=id_list[0])
        summary = Entrez.read(sum_handle)
        sum_handle.close()
        
        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        genomic_info = doc.get("GenomicInfo", [{}])[0]
        
        # Retorna o mapeamento preciso das coordenadas cromossômicas
        return {
            "gene": gene_symbol,
            "gene_id": id_list[0],
            "chr_acc": genomic_info.get("ChrAccVer", "N/A"),
            "chr_num": genomic_info.get("ChrLoc", "N/A"),
            "start": genomic_info.get("ChrStart", 0),
            "stop": genomic_info.get("ChrStop", 0)
        }

    def download_chromosome_fasta(self, accession_id: str, output_path: str) -> str:
        """
        Realiza o download streaming do arquivo FASTA do cromossomo via efetch REST.
        A escrita em blocos (chunks) evita o estouro de memória RAM com arquivos grandes.
        
        Parameters:
            accession_id (str): ID de acesso do cromossomo RefSeq (ex: 'NC_052532.1').
            output_path (str): Caminho onde o arquivo FASTA será salvo.
            
        Returns:
            str: Caminho do arquivo baixado.
        """
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nuccore",
            "id": accession_id,
            "rettype": "fasta",
            "retmode": "text"
        }
        print(f"\nIniciando download do FASTA ({accession_id})... Isso pode levar alguns minutos.")
        
        # Faz a requisição HTTP com stream ativado para leitura progressiva
        response = requests.get(url, params=params, stream=True)
        response.raise_for_status()
        
        # Grava o arquivo no disco em blocos de 1 MB
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                
        print(f"Download concluído: Salvo em {output_path}")
        return output_path

    @staticmethod
    def stream_fasta_efficiently(file_path: str):
        """
        Gerador (Generator) para leitura otimizada de arquivos FASTA pesados.
        Emite registro por registro via 'yield' para evitar o carregamento
        completo de cromossomos massivos na memória RAM.
        """
        print("Carregando arquivo FASTA via streaming (Memory-Friendly)...")
        for record in SeqIO.parse(file_path, "fasta"):
            yield record