# PKI: a bioinformatics method of quantifying the importance of nodes in gene regulatory network via a pseudo knockout index

## Highlights:
*	A method called Pseudo Knockout Importance (PKI) is proposed to quantify the importance of nodes and node combinations in gene regulatory network.
*	PKI designs gene pseudo knockout experiments and defines a criteria score based on the coefficient of determination.
*	PKI models dynamic gene regulations via ODE model and employs FPCA method to achieve an accurate estimation of gene expression differential value.
*	PKI formulates the key combinations with different numbers of nodes by a mathematical programming model and applies genetic algorithm to solve it efficiently.
*	The case studies in ESC differentiation generate promising alternative key TF combinations in reprogramming iPS cells.

## Workflow:
![image](https://user-images.githubusercontent.com/54654413/211350070-00aac877-6519-4cc3-82df-6b4fe3124ab3.png)
Figure 1. The workflow of PKI method. (A) ODE model describes cell-specific GRN. An ODE model is employed to describe the dynamic gene regulations over time. (B) Pseudo knockout rules. If node TF1 is supposed to be the knockout gene, edges linked to it and the isolated nodes will be removed simultaneously. (C) Pseudo knockout importance score. PKI score evaluation criteria is based on the coefficient of determination in regression, which reflects the consistency between gene expression data and ODE model. (D) Evaluation. The ranking of TFs can be sorted based on PKI scores. For TF combinations, genetic algorithm is applied to solve an optimization problem for obtaining gene combinations.

* If you have any questions about **PKI**, please directly contact the corresponding author [Prof. Zhi-Ping Liu](https://scholar.google.com/citations?user=zkBXb_kAAAAJ&hl=zh-CN&oi=ao) with the E-mail: zpliu@sdu.edu.cn

