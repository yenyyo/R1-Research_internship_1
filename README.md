\documentclass{article}
\usepackage{hyperref}
\usepackage{enumitem}

\title{Research Internship 1}
\author{Pablo de Vicente Abad}
\date{26th June 2024}

\begin{document}

\maketitle

\section*{Overview}
This repository contains the code and documentation created during my first research internship. The project focuses on immunosequencing to identify signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire, based on the work by Emerson et al.

\textbf{Note}: Folders preceded by '00 -' are auxiliary and not essential for understanding the main project.

\section*{Table of Contents}
\begin{enumerate}
    \item \hyperref[sec:installation]{Installation}
    \item \hyperref[sec:code_overview]{Code Overview}
    \begin{enumerate}[label*=\arabic*.]
        \item \hyperref[sec:data_preparation]{Data Preparation}
        \item \hyperref[sec:clustering]{Clustering}
        \item \hyperref[sec:analysis]{Analysis}
    \end{enumerate}
\end{enumerate}

\section*{Installation}
\label{sec:installation}
To get started, follow these steps:
\begin{enumerate}
    \item Clone the repository.
    \item Download the dataset used in the Emerson paper: \\
    \textit{Emerson, R., DeWitt, W., Vignali, M. et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat Genet 49, 659–665 (2017). \href{https://doi.org/10.1038/ng.3822}{https://doi.org/10.1038/ng.3822}}
    \item Follow the installation steps inside \texttt{00-Docu/01-Setup/}
\end{enumerate}

\section*{Code Overview}
\label{sec:code_overview}
\subsection*{Data Preparation}
\label{sec:data_preparation}
Preprocessing techniques applied to the original dataset can be found in \texttt{01-Code/data\_preparation/}:
\begin{enumerate}
    \item \texttt{cugraph\_clustering}: Original notebook.
    \item \texttt{data\_analysis}: Initial analysis of the dataset.
    \item \texttt{read\_emerson\_v0} / \texttt{aux\_funcs\_v0}: Generated the dataset used in later steps. TCRs were defined as a CDR3 sequence.
    \item \texttt{read\_emerson\_v2} / \texttt{aux\_funcs\_v2}: Generated the dataset used in later steps. TCRs were defined as a V\_gene + CDR3 sequence.
\end{enumerate}

\subsection*{Clustering}
\label{sec:clustering}
\begin{enumerate}
    \item \texttt{Create\_clusters}: Creates the clusters of TCRs based on co-occurrence.
    \item \texttt{Analyse\_clusters} / \texttt{aux\_clustering}: Loads data to compare and filter different clusters. Generates the file \texttt{gmm\_results.pkl}.
\end{enumerate}

\subsection*{Analysis}
\label{sec:analysis}
\begin{enumerate}
    \item \texttt{hla\_analysis}: Analyzes each of the groups that divide each cluster and extracts HLA differences between both.
\end{enumerate}

\end{document}
