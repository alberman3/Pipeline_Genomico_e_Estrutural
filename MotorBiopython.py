import re
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import Entrez, SeqIO

class MotorBiopython:
    def __init__(self, chr_acc, start, stop, fasta_file, gene_symbol="GAPDH"):
        self.chr_acc = chr_acc
        self.fasta_file = fasta_file
        self.gene_symbol = gene_symbol
        
        # Garante que o start é o menor número para a API de download
        self.genomic_start = min(int(start), int(stop)) + 1 
        self.genomic_stop = max(int(start), int(stop)) + 1
        
        self.exons_global = []
        self.is_reverse = False
        
    def executar_fase2(self):
        print("\n" + "="*50)
        print("--- INICIANDO FASE 2: MOTOR DE PROCESSAMENTO ---")
        print("="*50)
        
        # ==========================================
        # PASSO 3: Parse da Anotação e Plotagem
        # ==========================================
        print("\n1. Baixando Anotação GenBank do Locus e aplicando RegEx...")
        
        # Baixa apenas o trecho do cromossomo referente ao gene (Economiza memória)
        handle = Entrez.efetch(
            db="nuccore", id=self.chr_acc, rettype="gb", retmode="text", 
            seq_start=self.genomic_start, seq_stop=self.genomic_stop
        )
        gb_text = handle.read()
        handle.close()
        
        # Expressão Regular (re) para encontrar blocos "CDS" (Sequência Codificante)
        # Suporta tanto fita direta: join(1..10, 20..30) quanto fita reversa: complement(join(...))
        cds_match = re.search(r'CDS\s+(?:complement\()?join\((.*?)\)', gb_text, re.DOTALL)
        
        if not cds_match:
            # Tenta sem o 'join' (genes de éxon único)
            cds_match = re.search(r'CDS\s+(?:complement\()?(\d+\.\.\d+)\)?', gb_text)
            if not cds_match:
                 raise ValueError("Erro: Coordenadas CDS não encontradas na anotação.")
        
        # Descobre se está na fita negativa (3' -> 5') buscando a palavra 'complement'
        self.is_reverse = "complement" in cds_match.group(0)
        
        raw_coords = cds_match.group(1)
        # Extrai todos os pares de números (início..fim) ignorando quebras de linha
        exon_pairs = re.findall(r'(\d+)\.\.(\d+)', raw_coords)
        
        # Converte as coordenadas locais (do recorte) para as coordenadas globais do cromossomo
        for local_start, local_end in exon_pairs:
            global_s = self.genomic_start + int(local_start) - 1
            global_e = self.genomic_start + int(local_end) - 1
            self.exons_global.append((global_s, global_e))
            
        print(f"Éxons identificados no Cromossomo: {len(self.exons_global)}")
        print(f"Fita Direção: {'Reversa (3\' -> 5\')' if self.is_reverse else 'Direta (5\' -> 3\')'}")
        
        print("\n2. Gerando Gráfico (Figura 1)...")
        fig, ax = plt.subplots(figsize=(10, 2))
        ax.plot([self.genomic_start, self.genomic_stop], [0, 0], color="black", zorder=1, label="Gene Span (Introns)")
        
        for i, (ex_start, ex_end) in enumerate(self.exons_global):
            width = ex_end - ex_start
            rect = plt.Rectangle((ex_start, -0.2), width, 0.4, color="darkorange", zorder=2, label="Éxon (CDS)" if i==0 else "")
            ax.add_patch(rect)
            
        ax.set_ylim(-1, 1)
        ax.set_yticks([])
        ax.set_xlabel(f"Posição no Cromossomo {self.chr_acc} (bp)")
        ax.set_title(f"Distribuição Espacial de Éxons - {self.gene_symbol} (Gallus gallus)")
        ax.legend(loc="upper right")
        
        grafico_nome = f"Figura1_{self.gene_symbol}_exons.png"
        plt.tight_layout()
        plt.savefig(grafico_nome, dpi=300)
        print(f"Gráfico salvo com sucesso: '{grafico_nome}'")

        # ==========================================
        # PASSO 4: Slicing e Algoritmo de Fita Reversa
        # ==========================================
        print("\n3. Fazendo Slicing no arquivo FASTA massivo...")
        cds_sequence = ""
        
        # Abre o FASTA massivo e encontra o cromossomo sem estourar a memória RAM
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            if self.chr_acc in record.id:
                for ex_start, ex_end in self.exons_global:
                    # -1 porque o índice do Python começa em 0
                    cds_sequence += record.seq[ex_start - 1 : ex_end] 
                break
                
        dna_seq = Seq(str(cds_sequence))
        
        if self.is_reverse:
            print("Aplicando .reverse_complement() para corrigir a fita negativa...")
            dna_seq = dna_seq.reverse_complement()
            
        # ==========================================
        # PASSO 5: Tradução in silico
        # ==========================================
        print("\n4. Traduzindo DNA -> Proteína (.translate())...")
        protein_seq = dna_seq.translate()
        
        out_fasta = f"PROTEINA_{self.gene_symbol}_{self.taxon.replace(' ', '_')}.fasta"
        with open(out_fasta, "w") as f:
            f.write(f">{self.gene_symbol} | {self.taxon} | Extraido e Traduzido via Biopython\n")
            f.write(str(protein_seq) + "\n")
            
        print(f"Sucesso! Arquivo gerado: '{out_fasta}'")
        print(f"Tamanho da Proteína: {len(protein_seq)} aminoácidos.")
        print("="*50)