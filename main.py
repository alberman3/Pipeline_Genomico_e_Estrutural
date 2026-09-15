import os
from miner import NCBIGenomeMiner  # Importa a classe do arquivo miner.py
from MotorBiopython import MotorBiopython  # Importa a classe do motor da Fase 2

if __name__ == "__main__":
    # Instancia a ferramenta para a espécie designada
    miner = NCBIGenomeMiner(taxon="Gallus gallus")
    
    # -------------------------------------------------------------
    # PASSO 1: Mapeamento dos Cromossomos e geração da Tabela 1
    # -------------------------------------------------------------
    print("--- [PASSO 1] Obtendo metadados dos cromossomos ---")
    df_chromosomes = miner.get_chromosome_metadata()
    
    # Exibe a Tabela 1 solicitada no RTA
    print(df_chromosomes.head(10)) 
    df_chromosomes.to_csv("tabela1_caracterizacao_genoma_gallus_gallus.csv", index=False)
    
    # -------------------------------------------------------------
    # PASSO 2: Localização do Gene GAPDH e Download do Cromossomo
    # -------------------------------------------------------------
    print("\n--- [PASSO 2] Localizando o gene GAPDH ---")
    gene_info = miner.locate_gene(gene_symbol="GAPDH")
    print(f"Gene GAPDH localizado no Cromossomo RefSeq: {gene_info['chr_acc']} (Cromossomo nº {gene_info['chr_num']})")
    
    # Download do FASTA do cromossomo
    fasta_filename = f"chromosome_{gene_info['chr_acc']}.fasta"
    if not os.path.exists(fasta_filename):
        miner.download_chromosome_fasta(accession_id=gene_info["chr_acc"], output_path=fasta_filename)
    
    # -------------------------------------------------------------
    # FASE 2: Motor de Processamento Biopython (Passos 3, 4 e 5)
    # -------------------------------------------------------------
    # Processamento e leitura eficiente com Biopython (SeqIO)
    if os.path.exists(fasta_filename):
        motor = MotorBiopython(
            chr_acc=gene_info['chr_acc'],
            start=gene_info['start'],
            stop=gene_info['stop'],
            fasta_file=fasta_filename,
            taxon="Gallus gallus",
            gene_symbol="GAPDH"
        )
        motor.executar_fase2()