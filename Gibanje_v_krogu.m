function[p, pot, x, y] = Gibanje_v_krogu(R, smer, korak, ponovitve, l)
    %funkcija Gibanje_v_krogu s pomočjo statistične simulacije Monte Carlo 
    %izračuna verjetnost, da delec iz središča kroga pride skozi lok na
    %krožnici dolžine l. Delec se premika naključno (Brownovo gibanje). 
    
    %Parametri funkcije so:
    %R ... radij krožnice.
    %smer ... v koliko smeri se delec lahko giba. Določena je enakomerno na
    %intervalu od [0, 2π]. Če je smer npr. 4 se delec lahko giblje 
    % na π/2, π, 3π/2 ali 2π.
    %ponovitve ... število ponovitev opisanega gibanja. Ena ponovitev šteje, ko
    %je delec od središča oddaljen za Rd ≥ R.
    %l ... dolžina loka na krožnici.
    
    %Funkcija vrne količnik p med številom gibanj, ki so šli skozi lok l na
    %krožnici in številom vseh gibanj.
    
    interval = linspace(0, 360, smer + 1);

    smeri = interval(1:end - 1);
    
    zadetki = 0;
    x = 0;
    y = 0;

    for i = 1:ponovitve
        stevec = 2;
        pozicija = [0 0];
        pot = [0; 0]; %x in y
        while true
            %definirajmo smer
            smer = sym(deg2rad(smeri(randi(length(smeri)))));

            pozicija(1) = pozicija(1) + korak * cos(smer);
            pozicija(2) = pozicija(2) + korak * sin(smer);
            pot(1, end + 1) = pozicija(1);
            pot(2, stevec) = pozicija(2);
            if norm(pozicija) >= R
                %Daljica AB je sestavljena iz začetne točke A = (x1, y1) in
                %končne točke B = (x2, y2)
                x1 = pot(1, end - 1);
                y1 = pot(2, end - 1);
                x2 = pozicija(1);
                y2 = pozicija(2);

                %Iščemo koeficiente premice, ki seka krog.
              if x1 == x2
                  k = Inf;
                  n = x1;
              else
                  coefficients = polyfit([x1, x2], [y1, y2], 1);
                  k = coefficients (1);
                  n = coefficients (2);
              end
              %Poiščimo presečišče daljice AB s krogom.
              [P1,P2] = linecirc(k,n,0,0,R);
              %Preverimo katero presečišče je pravo.
              if (min(x1, x2) <= P1(1)) && (P1(1) <= max(x1, x2)) && (min(y1, y2) <= P2(1)) && (P2(1) <= max(y1, y2))
                  x = P1(1); 
                  y = P2(1);
              else
                  x = P1(2);
                  y = P2(2);
              end
              kot_P = mod(atan2(y,x),2*pi);

                %Preverimo ali je presečišče v našem želenem intervalu.
                if l >= kot_P && kot_P >= 0
                    zadetki = zadetki + 1;
                end
          
                break
            end
            stevec = stevec + 1;
        end
    end
    %izračunajmo še dobljeno verjetnost
    p = zadetki/ponovitve;
end
