for k=1:10:N-1
    plot(v(k+1,end),(v(k+1,end)-v(k,end))/(L/N),'bx')
    hold on
    pause
end