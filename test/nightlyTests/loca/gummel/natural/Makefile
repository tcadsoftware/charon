IVSWEEP=../../../utils/scripts/python/iv-sweep.py

default:
	@echo "You didn't speciy anything to do!"

clean:
	rm -f currents.dat currents-loca.dat pndiode.dd-result.[0-9]*.exo pndiode.nlp-result.exo pndiode*.txt iv-sweep-inp.xml 

run:
	$(MAKE) clean
	$(CHARONEXE) --i=pndiode.nlp.xml
	$(CHARONEXE) --current --i=pndiode.loca.dd.xml
	$(IVSWEEP) --vstart 0.0 --vinc 0.1 --nsteps 11 --charon ~/Projects/Charon2/builds/divcurl.gnu.opt/tcad-charon/charon_mp.exe --ig-filename pndiode.nlp-result.exo --save-screen-output --ncdump ~/Software/bin/ncdump pndiode.sweep.dd.xml.TMPL
