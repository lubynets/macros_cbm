#!/bin/bash

# iParticle
Particle=('Lambda' 'Kshort' 'Xi')
Particle_greek=('$\Lambda$' '$K^0_S$' '$\Xi^-$')
Particle_pdg=('3122' '310' '3312')

# iSetup
Setup=('d12' 'u12' 'd3')
Setup_title=('DCM-QGSM-SMM @ 12A GeV/c' 'UrQMD @ 12A GeV/c' 'DCM-QGSM-SMM @ 3.3A GeV/c' )
Evegen=('dcmqgsm' 'urqmd' 'dcmqgsm')
Beam_mom=('12' '12' '3.3')

# iResolution
RP_setup=('$\Psi_{RP}$'\
          'PSD1 with MC-true resolution' 'PSD2 with MC-true resolution' 'PSD3 with MC-true resolution'\
          'PSD1 with resolution determined with 3-subevents method' 'PSD2 with resolution determined with 3-subevents method' 'PSD3 with resolution determined with 3-subevents method'\
          'PSD1 with resolution determined with 4-subevents method' 'PSD2 with resolution determined with 4-subevents method' 'PSD3 with resolution determined with 4-subevents method')
RP_setup_short=('psirp'\
                'psd1mc' 'psd2mc' 'psd3mc'\
                'psd1sub3' 'psd2sub3' 'psd3sub3'\
                'psd1sub4' 'psd2sub4' 'psd3sub4')

for iSetup in {0..2} # Setup (d12 u12 d3)
do
for iResolution in {0..9} # RP_setup
do
for iParticle in {0..2} # Particle
do
SLOPE_IN=${Setup[iSetup]}/${Particle[iParticle]}/dv1dy/slope.${Evegen[iSetup]}.${Beam_mom[iSetup]}agev.${Particle_pdg[iParticle]}.pdf
INTERCEPT_IN=${Setup[iSetup]}/${Particle[iParticle]}/dv1dy/intercept.${Evegen[iSetup]}.${Beam_mom[iSetup]}agev.${Particle_pdg[iParticle]}.pdf
ls -l $SLOPE_IN
ls -l $INTERCEPT_IN

SLOPE_OUT=dv1dy.slope.${Setup[iSetup]}.${Particle[iParticle]}.${RP_setup_short[iResolution]}.pdf
INTERCEPT_OUT=dv1dy.intercept.${Setup[iSetup]}.${Particle[iParticle]}.${RP_setup_short[iResolution]}.pdf

# pdftk $SLOPE_IN cat $((iResolution+1)) output $SLOPE_OUT
# pdftk $INTERCEPT_IN cat $((iResolution+1)) output $INTERCEPT_OUT


if [[ $(($iParticle % 3)) -eq 0 ]]
then
echo \\begin{figure}[h!] >> file.tex
echo \\centering >> file.tex
fi
echo \\includegraphics[width=0.46\\textwidth]{$SLOPE_OUT} >> file.tex
echo \\hspace{0.04\\textwidth} >> file.tex
echo \\includegraphics[width=0.46\\textwidth]{$INTERCEPT_OUT} >> file.tex
if [[ $(($iParticle % 3)) -eq 2 ]]
then
echo \\\\ >> file.tex
echo \\caption{Centrality dependence of \$dv_1/dy\$ \(\\textbf\{left\}\) slope and \(\\textbf\{right\}\) intercept of \(\\textbf\{top\}\) ${Particle_greek[0]}, \(\\textbf\{middle\}\) ${Particle_greek[1]}, \(\\textbf\{bottom\}\) ${Particle_greek[2]}. >> file.tex
echo Input generated with ${Setup_title[iSetup]}, reaction plane determined from ${RP_setup[iResolution]}.} >> file.tex
echo \\label{fig:dv1dy:${Setup[iSetup]}:${RP_setup_short[iResolution]}} >> file.tex
echo \\end{figure} >> file.tex
echo \\newpage >> file.tex
echo >> file.tex
fi

done
done
done
