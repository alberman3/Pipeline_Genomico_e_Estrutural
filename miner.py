import os
import requests
import pandas as pd
from Bio import Entrez, SeqIO
from typing import Dict

class NCBIGenomeMiner:
    """
    Minerador Genômico Refatorado para usar Bio.Entrez (Via Oficial e Segura do NCBI).
    """
    def __init__(self, taxon: str):
        self.taxon = taxon
        
        # Obrigatório: O NCBI exige que você se identifique para não bloquear seu IP
        Entrez.email = "albermangomes3@gmail.com"
        Entrez.tool = "ProjetoBioinfo2026"

    def get_chromosome_metadata(self) -> pd.DataFrame:
        print(f"--- [PASSO 1] Consultando genoma de {self.taxon} via Entrez ---")
        chromosomes_data = []

        try:
            # 1. Procura a montagem (Assembly) de referência da espécie
            print("1. Buscando o genoma de referência (Assembly)...")
            # Etapa 1: Busca o genoma de referência da espécie
            search_handle = Entrez.esearch(db="assembly", term=f'"{self.taxon}"[Organism] AND "latest refseq"[filter]', retmax=1)
            
            assembly_record = Entrez.read(search_handle)
            search_handle.close()

            if not assembly_record.get("IdList"):
                raise ValueError("Nenhum genoma encontrado. Verifique o nome da espécie.")

            # 2. Obtém o número de acesso da montagem (ex: GCF_016699485.2)
            sum_handle = Entrez.esummary(db="assembly", id=assembly_record["IdList"][0])
            assembly_summary = Entrez.read(sum_handle)
            sum_handle.close()
            
            assembly_acc = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            print(f"Montagem encontrada com sucesso: {assembly_acc}")

            # 3. Mapeia os cromossomos ligados a esta montagem
            print("2. Mapeando os cromossomos (Nuccore)...")
            # Etapa 2: Mapeia os cromossomos
            nuc_search = Entrez.esearch(db="nuccore", term=f'{assembly_acc}[Assembly] AND biomol_genomic[PROP]', retmax=1500)
            nuc_record = Entrez.read(nuc_search)
            nuc_search.close()
            
            seq_ids = nuc_record.get("IdList", [])
            
            # 4. Coleta o tamanho de cada sequência
            nuc_sum = Entrez.esummary(db="nuccore", id=",".join(seq_ids))
            seq_summaries = Entrez.read(nuc_sum)
            nuc_sum.close()

            for seq in seq_summaries:
                acc = seq.get("AccessionVersion", "")
                length = seq.get("Length", 0)
                title = seq.get("Title", "")
                
                if length and int(length) > 1:
                    # Verifica no título INTEIRO se é um cromossomo
                    if "chromosome" in title.lower():
                        # Tenta extrair o nome limpo (ex: "chromosome 1")
                        # Procura algo como "chromosome X" no texto
                        import re
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

        df = pd.DataFrame(chromosomes_data)
        if df.empty:
            raise RuntimeError(f"Falha definitiva ao extrair cromossomos para {self.taxon}.")
            
        df = df.drop_duplicates(subset=["ID do Cromossomo (FASTA)"])
        # Retorna ordenado do maior para o menor
        return df.sort_values(by="Tamanho (bp)", ascending=False).reset_index(drop=True)

    def locate_gene(self, gene_symbol: str) -> Dict[str, str]:
        print(f"\n--- [PASSO 2] Localizando coordenadas do gene {gene_symbol} ---")
        #Aqui nós precisamos descobrir exatamente em qual dos 39 pares de cromossomos da galinha o gene GAPDH está escondido.
        search_handle = Entrez.esearch(db="gene", term=f'"{self.taxon}"[Organism] AND {gene_symbol}[Gene Name]', retmax=1)
        record = Entrez.read(search_handle)
        search_handle.close()
        
        id_list = record.get("IdList", [])
        if not id_list:
            raise ValueError(f"Gene '{gene_symbol}' não encontrado para {self.taxon}.")
            
        sum_handle = Entrez.esummary(db="gene", id=id_list[0])
        summary = Entrez.read(sum_handle)
        sum_handle.close()
        
        doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
        genomic_info = doc.get("GenomicInfo", [{}])[0]
        
        return {
            "gene": gene_symbol,
            "gene_id": id_list[0],
            "chr_acc": genomic_info.get("ChrAccVer", "N/A"),
            "chr_num": genomic_info.get("ChrLoc", "N/A"),
            "start": genomic_info.get("ChrStart", 0),
            "stop": genomic_info.get("ChrStop", 0)
        }

    def download_chromosome_fasta(self, accession_id: str, output_path: str) -> str:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nuccore",
            "id": accession_id,
            "rettype": "fasta",
            "retmode": "text"
        }
        print(f"\nIniciando download do FASTA ({accession_id})... Isso pode levar alguns minutos.")
        
        response = requests.get(url, params=params, stream=True)
        response.raise_for_status()
        
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                
        print(f"Download concluído: Salvo em {output_path}")
        return output_path

    @staticmethod
    def stream_fasta_efficiently(file_path: str):
        print("Carregando arquivo FASTA via streaming (Memory-Friendly)...")
        for record in SeqIO.parse(file_path, "fasta"):
            yield record