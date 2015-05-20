%%% ====================================================================
%%% @LaTeX-file{
%%%    filename  = "mcom-l.cls",
%%%    version   = "2.01",
%%%    date      = "2007/06/18",
%%%    time      = "09:17:48 EDT",
%%%    checksum  = "31108 41 149 1479",
%%%    author    = "Michael Downes",
%%%    address   = "American Mathematical Society, PO Box 6248,
%%%                       Providence, RI 02940, USA",
%%%    email     = "tech-support@ams.org (Internet)",
%%%    supported = "yes",
%%%    abstract  = "This is a LaTeX documentclass for the AMS journal
%%%                 Mathematics of Computation.",
%%%    docstring = "The checksum field above contains a CRC-16 checksum
%%%                 as the first value, followed by the equivalent of
%%%                 the standard UNIX wc (word count) utility output
%%%                 of lines, words, and characters.  This is produced
%%%                 by Robert Solovay's checksum utility.",
%%% }
%%% ====================================================================

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mcom-l}
              [2007/06/18 v2.01 MCOM Author Class]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{amsart}}
\ProcessOptions\relax

\LoadClass{amsart}[1996/10/24]

\def\publname{MATHEMATICS OF COMPUTATION}

\def\ISSN{0025-5718}

\copyrightinfo{\currentyear}%    % copyright year
  {American Mathematical Society}% copyright holder

\renewcommand{\datename}{Received by the editor}

\endinput
