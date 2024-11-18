package an;

import java.util.Scanner;

public class Semestral {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        System.out.print("Digite a ordem da matriz (n): ");
        int n = scanner.nextInt();
        double[][] A = new double[n][n];
        double[] B = new double[n];

        System.out.println("Digite a matriz A (coeficientes):");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                System.out.print("A[" + i + "][" + j + "]: ");
                A[i][j] = scanner.nextDouble();
            }
        }

        System.out.println("Digite o vetor B (termos independentes):");
        for (int i = 0; i < n; i++) {
            System.out.print("B[" + i + "]: ");
            B[i] = scanner.nextDouble();
        }

        // Transformar o sistema na forma x = Gx + H
        double[][] G = new double[n][n];
        double[] H = new double[n];
        transformarSistema(A, B, G, H);

        // Calcular a norma de G
        double normaG = calcularNormaMatriz(G);
        System.out.printf("Norma de G: %.5f\n", normaG);

        if (normaG < 1) {
            System.out.print("Digite a precisão desejada (epsilon): ");
            double epsilon = scanner.nextDouble();

            System.out.print("Digite o vetor inicial X0 (separado por espaços): ");
            double[] x0 = new double[n];
            for (int i = 0; i < n; i++) {
                x0[i] = scanner.nextDouble();
            }

            // Resolver o sistema com o vetor inicial fornecido e imprimir a tabela de iterações
            double[][] iteracoes = resolverSistema(G, H, epsilon, x0);

            if (iteracoes != null) {
                String format = "%." + casasDecimais(epsilon) + "f";
                System.out.println("\nTabela de iterações:");
                System.out.printf("%-5s %-10s %-10s %-10s %-20s %-10s\n", "Iter", "x1", "x2", "x3", "epsilon||Xm+1-Xm||", "Status");
                
                System.out.printf("%-5d " + format + " " + format + " " + format + " %-20s %-10s\n", 0, x0[0], x0[1], x0[2], " - - - - -", "Não");
                
                for (int i = 0; i < iteracoes.length; i++) {
                    System.out.printf("%-5d ", i + 1); // Iterações começam em 1
                    for (int j = 0; j < n; j++) {
                        System.out.printf(format + " ", iteracoes[i][j]);
                    }
                    System.out.printf("%-20.5f ", iteracoes[i][n]);
                    System.out.printf("%-10s\n", iteracoes[i][n + 1] == 1.0 ? "Sim" : "Não");
                }

                // Mostrar o resultado final
                double x1 = iteracoes[iteracoes.length - 1][0];
                double x2 = iteracoes[iteracoes.length - 1][1];
                double x3 = iteracoes[iteracoes.length - 1][2];
                System.out.printf("\nResultado: x1=" + format + " x2=" + format + " x3=" + format + "\n", x1, x2, x3);
                System.out.printf("Solucao:\n" + format + " +/- " + format + "\n" + format + " +/- " + format + "\n" + format + " +/- " + format + "\n", x1, epsilon, x2, epsilon, x3, epsilon);
            } else {
                System.out.println("O sistema não possui solução ou não convergiu.");
            }
        } else {
            System.out.println("A norma de G é maior ou igual a 1. O método de Jacobi pode não convergir.");
        }

        scanner.close();
    }

    // Método para transformar o sistema na forma x = Gx + H
    public static void transformarSistema(double[][] A, double[] B, double[][] G, double[] H) {
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                if (i != j) {
                    G[i][j] = -A[i][j] / A[i][i];
                } else {
                    G[i][j] = 0;
                }
            }
            H[i] = B[i] / A[i][i];
        }

        // Exibir o sistema na forma x = Gx + H
        System.out.println("\nSistema de equações na forma x = Gx + H:");
        for (int i = 0; i < G.length; i++) {
            System.out.print("x" + (i + 1) + " = ");
            for (int j = 0; j < G[i].length; j++) {
                if (G[i][j] != 0) {
                    System.out.printf("%.5f*x%d ", G[i][j], (j + 1));
                }
            }
            System.out.printf("+ %.5f\n", H[i]);
        }
    }

    // Método para resolver o sistema de equações usando o Método de Jacobi
    public static double[][] resolverSistema(double[][] G, double[] H, double epsilon, double[] x0) {
        int n = H.length;
        double[] xAtual = x0.clone(); // Utiliza o vetor inicial fornecido
        double[] xNovo = new double[n];
        double[][] iteracoes = new double[1000][n + 2]; // Limite de 1000 iterações
        int iteracao = 0;
        boolean convergiu = false;

        while (!convergiu && iteracao < 1000) {
            for (int i = 0; i < n; i++) {
                double soma = 0;
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        soma += G[i][j] * xAtual[j];
                    }
                }
                xNovo[i] = soma + H[i];
            }

            double norma = calcularNormaVetor(xNovo, xAtual);
            System.arraycopy(xNovo, 0, iteracoes[iteracao], 0, n);
            iteracoes[iteracao][n] = norma;
            iteracoes[iteracao][n + 1] = (norma < epsilon) ? 1.0 : 0.0; // Status

            if (norma < epsilon) {
                convergiu = true;
            }

            System.arraycopy(xNovo, 0, xAtual, 0, n);
            iteracao++;
        }

        if (!convergiu) {
            System.out.println("Método não convergiu após 1000 iterações.");
            return null;
        }

        double[][] iteracoesFinais = new double[iteracao][n + 2];
        System.arraycopy(iteracoes, 0, iteracoesFinais, 0, iteracao);
        return iteracoesFinais;
    }

    // Função para calcular a norma ||Xm+1 - Xm|| para vetores
    public static double calcularNormaVetor(double[] xNovo, double[] xAtual) {
        double soma = 0;
        for (int i = 0; i < xNovo.length; i++) {
            soma += Math.pow(xNovo[i] - xAtual[i], 2);
        }
        return Math.sqrt(soma);
    }

    // Função para calcular a norma de uma matriz usando ||G|| = sqrt(somatório(G[i,j]^2))
    public static double calcularNormaMatriz(double[][] G) {
        double soma = 0;
        System.out.println("Detalhes do cálculo da norma de G:");
        for (double[] G1 : G) {
            for (int j = 0; j < G1.length; j++) {
                soma += G1[j] * G1[j];
            }
        }
        return Math.sqrt(soma);
    }

    // Função para contar as casas decimais de um número
    public static int casasDecimais(double numero) {
        String texto = Double.toString(Math.abs(numero));
        int indice = texto.indexOf(".");
        return texto.length() - indice - 1;
    }
}
