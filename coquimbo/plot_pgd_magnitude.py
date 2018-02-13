
plot(tpgd,Mw_nov)
plot(tpgd,ones(len(tpgd))*8.3,'--',c='r',lw=2)
ylim([6.5,9])
xlim([0,200])
xlabel('Seconds after OT')
ylabel('Magnitude')
grid()
errorbar(tpgd,Mw_nov,yerr=0.27,fmt='o',color='#6495ED',ecolor='#6495ED',markersize=3)