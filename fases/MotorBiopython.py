import os
import re
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import Entrez, SeqIO

class MotorBiopython:
    def __init__(self, chr_acc, start, stop, taxon, fasta_file, gene_symbol="GAPDH", output_dir="resultados/fase2"):
        self.chr_acc = chr_acc
        self.fasta_file = fasta_file
        self.gene_symbol = gene_symbol
        self.taxon = taxon
        self.output_dir = output_dir
        
        # Garante que a pasta de destino exista
        os.makedirs(self.output_dir, exist_ok=True)

        # Garante que o start é o menor número
        self.genomic_start = min(int(start), int(stop)) + 1 
        self.genomic_stop = max(int(start), int(stop)) + 1
        
        self.exons_global = []
        self.is_reverse = False

    def verificar_pertencimento_otimizado(self, cds_seq_str):
        """
        [OPÇÃO 1 OTIMIZADA] Validação instantânea via busca exata sem estourar a memória RAM.
        Examina apenas o trecho genômico e verifica os éxons individualmente.
        """
        print("\n--- [VALIDAÇÃO OTIMIZADA] Verificando pertencimento ao Cromossomo ---")
        
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            if self.chr_acc in record.id:
                # Extrai apenas o trecho do cromossomo referente ao locus do gene
                trecho_cromossomo = str(record.seq[self.genomic_start - 1 : self.genomic_stop]).upper()
                
                # Valida se os éxons individuais estão contidos no trecho recortado
                exons_encontrados = 0
                for ex_start, ex_end in self.exons_global:
                    # Ajusta coordenadas para o índice do recorte local
                    l_start = ex_start - self.genomic_start
                    l_end = ex_end - self.genomic_start
                    exon_seq = str(record.seq[ex_start - 1 : ex_end]).upper()
                    
                    if exon_seq in trecho_cromossomo or str(Seq(exon_seq).reverse_complement()) in trecho_cromossomo:
                        exons_encontrados += 1
                        
                if exons_encontrados == len(self.exons_global):
                    print(f"✓ CONFIRMADO: Todos os {exons_encontrados} éxons foram validados com 100% de precisão no Cromossomo!")
                    return True
                else:
                    print(f"⚠ AVISO: Apenas {exons_encontrados}/{len(self.exons_global)} éxons foram confirmados.")
                    return False
                    
        print("✗ ERRO: ID do Cromossomo não encontrado no arquivo FASTA.")
        return False

    def executar_fase2(self):
        print("\n" + "="*50)
        print("--- INICIANDO FASE 2: MOTOR DE PROCESSAMENTO ---")
        print("="*50)
        
        # ==========================================
        # PASSO 3: Parse da Anotação e Plotagem
        # ==========================================
        print("\n1. Baixando Anotação GenBank do Locus e aplicando RegEx...")
        
        handle = Entrez.efetch(
            db="nuccore", id=self.chr_acc, rettype="gb", retmode="text", 
            seq_start=self.genomic_start, seq_stop=self.genomic_stop
        )
        gb_text = handle.read()
        handle.close()
        
        cds_match = re.search(r'CDS\s+(?:complement\()?join\((.*?)\)', gb_text, re.DOTALL)
        if not cds_match:
            cds_match = re.search(r'CDS\s+(?:complement\()?(\d+\.\.\d+)\)?', gb_text)
            if not cds_match:
                 raise ValueError("Erro: Coordenadas CDS não encontradas na anotação.")
        
        self.is_reverse = "complement" in cds_match.group(0)
        raw_coords = cds_match.group(1)
        exon_pairs = re.findall(r'(\d+)\.\.(\d+)', raw_coords)
        
        for local_start, local_end in exon_pairs:
            global_s = self.genomic_start + int(local_start) - 1
            global_e = self.genomic_start + int(local_end) - 1
            self.exons_global.append((global_s, global_e))
            
        print(f"Éxons identificados no Cromossomo: {len(self.exons_global)}")
        print(f"Fita Direção: {'Reversa (3\' -> 5\')' if self.is_reverse else 'Direta (5\' -> 3\')'}")
        
        print("\n2. Gerando Gráfico (Figura 1)...")
        fig, ax = plt.subplots(figsize=(12, 3.5))
        
        # Extensão total do locus (de 0 bp ao tamanho final do gene)
        tamanho_total = self.genomic_stop - self.genomic_start
        
        # Desenha a linha dos íntrons iniciando em 0 bp
        ax.plot([0, tamanho_total], [0, 0], color="black", zorder=1, label="Gene Span (Introns)")
        
        # Estrutura para salvar o relatório de éxons
        relatorio_exons = []
        relatorio_exons.append(f"=== MAPEAMENTO DETALHADO DOS ÉXONS: {self.gene_symbol} ({self.taxon}) ===")
        relatorio_exons.append(f"Cromossomo: {self.chr_acc}")
        relatorio_exons.append(f"Orientação da Fita: {'Reversa (3\' -> 5\')' if self.is_reverse else 'Direta (5\' -> 3\')'}\n")
        relatorio_exons.append(f"{'Éxon':<8} | {'Início (bp)':<12} | {'Fim (bp)':<12} | {'Tamanho (bp)':<12}")
        relatorio_exons.append("-" * 52)

        for i, (ex_start, ex_end) in enumerate(self.exons_global):
            width = ex_end - ex_start
            
            # Posição relativa ao início do gene (Base 0)
            rel_start = ex_start - self.genomic_start
            
            # Desenha o bloco do Éxon no gráfico usando coordenadas relativas
            rect = plt.Rectangle((rel_start, -0.2), width, 0.4, color="darkorange", zorder=2, label="Éxon (CDS)" if i==0 else "")
            ax.add_patch(rect)
            
            # Texto acima do bloco: Identificação (E1, E2...) + Quantidade de pares de bases (Xbp)
            centro_x = rel_start + (width / 2)
            ax.text(centro_x, 0.3, f"E{i+1}\n({width}bp)", ha='center', va='bottom', fontsize=8, fontweight='bold', color='darkblue')
            
            # Adiciona ao relatório de texto (mantém as coordenadas cromossômicas globais)
            relatorio_exons.append(f"E{i+1:<7} | {ex_start:<12} | {ex_end:<12} | {width} bp")

        # ----------------------------------------------------
        # CÁLCULO E SOMA DOS ÉXONS (COM E SEM E1/E11)
        # ----------------------------------------------------
        tamanhos = [end - start for start, end in self.exons_global]
        soma_todos = sum(tamanhos)
        div_todos = soma_todos / 3

        tamanhos_filtrados = tamanhos[1:-1]
        soma_filtrados = sum(tamanhos_filtrados)
        div_filtrados = soma_filtrados / 3

        relatorio_exons.append("-" * 52)
        relatorio_exons.append(f"Soma de TODOS os Éxons (E1 a E{len(self.exons_global)}): {soma_todos} bp | Divisão por 3: {div_todos:.2f}")
        relatorio_exons.append(f"Soma sem 1º e Último (E2 a E{len(self.exons_global)-1}): {soma_filtrados} bp | Divisão por 3: {div_filtrados:.2f}")

        # Configuração do gráfico com escala relativa limpa
        plt.xticks(fontsize=8)
        
        ax.set_xlim(-100, tamanho_total + 100)
        ax.set_ylim(-0.8, 1.5)
        ax.set_yticks([])
        ax.set_xlabel(f"Posição Relativa no Gene (bp) [Locus Genômico: {self.chr_acc}:{self.genomic_start}-{self.genomic_stop}]")
        ax.set_title(f"Distribuição Espacial de Éxons - {self.gene_symbol} ({self.taxon})")
        ax.legend(loc="upper right")
        
        grafico_nome = os.path.join(self.output_dir, f"Figura1_{self.gene_symbol}_exons.png")
        plt.tight_layout()
        plt.savefig(grafico_nome, dpi=300)
        plt.close()
        print(f"Gráfico salvo com sucesso: '{grafico_nome}'")

        # ----------------------------------------------------
        # PRINT NO TERMINAL E SALVAMENTO EM ARQUIVO TXT
        # ----------------------------------------------------
        texto_completo = "\n".join(relatorio_exons)
        print("\n" + texto_completo + "\n")
        
        txt_nome = os.path.join(self.output_dir, f"coordenadas_exons_{self.gene_symbol}.txt")
        with open(txt_nome, "w") as f:
            f.write(texto_completo)
        print(f"Coordenadas dos éxons salvas no arquivo: '{txt_nome}'")

        # ==========================================
        # PASSO 4: Slicing e Algoritmo de Fita Reversa
        # ==========================================
        print("\n3. Fazendo Slicing no arquivo FASTA massivo...")
        cds_sequence = ""
        
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            if self.chr_acc in record.id:
                for ex_start, ex_end in self.exons_global:
                    cds_sequence += record.seq[ex_start - 1 : ex_end] 
                break
                
        dna_seq = Seq(str(cds_sequence))
        
        if self.is_reverse:
            print("Aplicando .reverse_complement() para corrigir a fita negativa...")
            dna_seq = dna_seq.reverse_complement()

        # ----------------------------------------------------
        # VALIDAÇÃO LEVE E INSTANTÂNEA
        # ----------------------------------------------------
        self.verificar_pertencimento_otimizado(str(dna_seq))

        # ----------------------------------------------------
        # SALVANDO A SEQUÊNCIA DE DNA (CDS) ANTES DA TRADUÇÃO
        # ----------------------------------------------------
        out_dna_fasta = os.path.join(self.output_dir, f"CDS_{self.gene_symbol}_{self.taxon.replace(' ', '_')}.fasta")
        with open(out_dna_fasta, "w") as f:
            f.write(f">{self.gene_symbol}_CDS | {self.taxon} | Sequencia Codificante (DNA 5'->3')\n")
            f.write(str(dna_seq) + "\n")
            
        print(f"Sucesso! Sequência de DNA (CDS) salva em: '{out_dna_fasta}'")
        print(f"Tamanho do CDS: {len(dna_seq)} bp.")

        # ----------------------------------------------------
        #  EXTRAIR EXONS (SEM O 1º E O ÚLTIMO)
        # ----------------------------------------------------
        print("\n[TAREFA EXTRA] Extraindo FASTA dos éxons isolados (limpando o 1º e o último)...")
        
        # Filtra a lista removendo o primeiro (0) e o último (-1)
        exons_filtrados = self.exons_global[1:-1]
        
        # Se for fita reversa, a ordem biológica dos éxons é lida de trás pra frente
        if self.is_reverse:
            exons_filtrados = list(reversed(exons_filtrados))
            
        out_exons_fasta = os.path.join(self.output_dir, f"EXONS_FILTRADOS_{self.gene_symbol}_{self.taxon.replace(' ', '_')}.fasta")
        
        with open(out_exons_fasta, "w") as f_out:
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                if self.chr_acc in record.id:
                    # Percorre apenas os éxons intermediários (do 2 até o penúltimo)
                    for i, (ex_start, ex_end) in enumerate(exons_filtrados):
                        # Índice ajustado para o cabeçalho (como pulamos o 1º, começamos do 2)
                        num_exon = i + 2 
                        
                        # Extrai a sequência exata do éxon no cromossomo
                        exon_seq = record.seq[ex_start - 1 : ex_end]
                        
                        # Se for fita negativa, aplica o complemento reverso
                        if self.is_reverse:
                            exon_seq = exon_seq.reverse_complement()
                            
                        # Escreve no formato Multi-FASTA pedido no quadro
                        f_out.write(f">Exon_{num_exon} | {self.gene_symbol} | {self.taxon}\n")
                        f_out.write(f"{str(exon_seq)}\n")
                    break
                    
        print(f"Sucesso! Éxons filtrados salvos no arquivo: '{out_exons_fasta}'")
            
        # ==========================================
        # PASSO 5: Tradução in silico
        # ==========================================
        print("\n4. Traduzindo DNA -> Proteína (.translate())...")
        protein_seq = dna_seq.translate()
        
        out_fasta = os.path.join(self.output_dir, f"PROTEINA_{self.gene_symbol}_{self.taxon.replace(' ', '_')}.fasta")
        with open(out_fasta, "w") as f:
            f.write(f">{self.gene_symbol} | {self.taxon} | Extraido e Traduzido via Biopython\n")
            f.write(str(protein_seq) + "\n")
            
        print(f"Sucesso! Arquivo gerado: '{out_fasta}'")
        print(f"Tamanho da Proteína: {len(protein_seq)} aminoácidos.")
        print("="*50)