#!/bin/bash

Particle=('lambda' 'kshort' 'xi' 'pipos' 'pineg')
Particle_greek=('$\Lambda$' '$K^0_S$' '$\Xi^-$' '$\pi^+$' '$\pi^-$')
Subevent=('psd1' 'psd2' 'psd3'\
          'etacut_1_charged' 'etacut_2_charged' 'etacut_3_charged'\
          'etacut_1_all' 'etacut_2_all' 'etacut_3_all')

for I in {0..4}
do
for J in {0..8}
do
# ls -l dv1dy.${Particle[I]}.${Subevent[J]}.pdf
# pdftk dv1dy.${Particle[I]}.${Subevent[J]}.pdf cat 1 output dv1dy.slope.${Particle[I]}.${Subevent[J]}.pdf
# pdftk dv1dy.${Particle[I]}.${Subevent[J]}.pdf cat 2 output dv1dy.offset.${Particle[I]}.${Subevent[J]}.pdf

if [[ $(($J % 3)) -eq 0 ]]
then
echo \\begin{figure}[h!] >> file.tex
fi
echo \\includegraphics[width=0.48\\textwidth]{dv1dy.slope.${Particle[I]}.${Subevent[J]}.pdf} >> file.tex
echo \\hspace{0.02\\textwidth} >> file.tex
echo \\includegraphics[width=0.48\\textwidth]{dv1dy.offset.${Particle[I]}.${Subevent[J]}.pdf} >> file.tex
if [[ $(($J % 3)) -eq 2 ]]
then
echo \\caption{\$dv_1/dy\$ \(\\textbf\{left\}\) slope and \(\\textbf\{right\}\) intercept of ${Particle_greek[I]}} >> file.tex
echo \\label{fig:${Particle[I]}} >> file.tex
echo \\end{figure} >> file.tex
echo \\newpage >> file.tex
echo >> file.tex
fi
done
done
