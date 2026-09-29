import os
from fases.miner import NCBIGenomeMiner  # Importa a classe do pasta fases/
from fases.MotorBiopython import MotorBiopython  # Importa a classe do motor da Fase 2
from fases.fase3_integracao import GAPDHFunctionalAPI

if __name__ == "__main__":
    # Garante a existência das pastas de resultados antes de salvar
    os.makedirs("resultados/fase1", exist_ok=True)
    os.makedirs("resultados/fase2", exist_ok=True)
    os.makedirs("resultados/fase3", exist_ok=True)

    # Instancia a ferramenta para a espécie designada
    miner = NCBIGenomeMiner(taxon="Gallus gallus")
    
    # -------------------------------------------------------------
    # PASSO 1: Mapeamento dos Cromossomos e geração da Tabela 1
    # -------------------------------------------------------------
    print("--- [PASSO 1] Obtendo metadados dos cromossomos ---")
    df_chromosomes = miner.get_chromosome_metadata()
    
    # Exibe a Tabela 1 solicitada no RTA
    print(df_chromosomes.head(10)) 
    tabela1_path = "resultados/fase1/tabela1_caracterizacao_genoma_gallus_gallus.csv"
    df_chromosomes.to_csv(tabela1_path, index=False)
    print(f"[OK] Tabela 1 salva em '{tabela1_path}'.")
    
    # -------------------------------------------------------------
    # PASSO 2: Localização do Gene GAPDH e Download do Cromossomo
    # -------------------------------------------------------------
    print("\n--- [PASSO 2] Localizando o gene GAPDH ---")
    gene_info = miner.locate_gene(gene_symbol="GAPDH")
    print(f"Gene GAPDH localizado no Cromossomo RefSeq: {gene_info['chr_acc']} (Cromossomo nº {gene_info['chr_num']})")
    
    # Download do FASTA do cromossomo para a pasta da fase1
    fasta_filename = f"resultados/fase1/chromosome_{gene_info['chr_acc']}.fasta"
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

    # -------------------------------------------------------------
    # FASE 3: Integração de APIs REST Funcionais (Passo 7)
    # -------------------------------------------------------------
    print("\n--- [PASSO 7 - FASE 3] Consultando APIs REST (UniProt e KEGG) ---")
    api_runner = GAPDHFunctionalAPI(uniprot_id="P00356")
    df_tabela2 = api_runner.generate_tabela_2()
    
    print("\n=== TABELA 2: DADOS FUNCIONAIS EXTRAÍDOS POR API (RTA 4.1) ===")
    print(df_tabela2.to_string(index=False))
    
    tabela2_path = "resultados/fase3/tabela2_dados_funcionais_gallus_gallus.csv"
    df_tabela2.to_csv(tabela2_path, index=False)
    print(f"\n[OK] Tabela 2 guardada em '{tabela2_path}'.")
    print("\n=== PIPELINE INDIVIDUAL CONCLUÍDO COM SUCESSO! ===")