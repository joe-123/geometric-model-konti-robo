function[robotShape] = GeometricModel(TACR,q)

% input: struct TACR, q:[jx3] configurational parameters
%
% output: robotShape.diskPoints: n rows for n disks, 12 columns for coordinates
%        (x,y,z) for central backbone, points tendon 1 (x,y,z), points tendon 2 (x,y,z),
%         points tendon 3 (x,y,z)
%         robotShape.diskRotation: n rows for n disks, columns represent
%         rotation matrices (3x3): columns 1-3 are matrix elements (1,1; 1,2; 1,3),
%         columns 4-6 are matrix elements (2,1; 2,2; 2,3), 
%         columns 7-9 are matrix elements (3,1; 3,2; 3,3),


% Die vorligende Funktion ist für eine beliebige Anzahl an Segmenten
% lauffaehig. Sonderfaelle treten nur beim ersten Segment auf, da dieses
% keine Anbindung aufweist. 
% Gruppe: Anton Sauer, Johannes Schulz
%%  
    % Ueberpruefung der Eingabe
    if abs(sum(q(:))) > 0.000001 % da double-Werte wird mit Minimum überprüft
        disp('Achtung: Summe der Seillaengen ist nicht Null!');
    end
    
    % Auswertung der uebergebenen Parameter
    number_of_segments = size(q, 1);
    number_of_discs = sum(TACR.ndisks);
    
    % Init der Variablen mit Startwerten
    diskHeightEnd = 0;
    theta_old = pi/2;
    speicher = 1;
    number_of_tendons = 3;
    phi_total = 0;
    beta = 2*pi/number_of_tendons;
    rot_old = eye(3);
    p_tendon= zeros(3,3);
    points_old = zeros(3,1);
    robotShape.diskPoints = zeros(number_of_discs, 12);
    robotShape.diskRotation = zeros(number_of_discs, 9);
    
    % Schleife Ueber die Segmente
    for j=1:number_of_segments
        phi = -atan2(-q(j,1)*cos(beta)+q(j,2), -q(j,1)*sin(beta)); % Da phi negativ laeuft wird Vorzeichen umgedreht
        phi_total = phi_total+phi; % phi wird aufsummiert um delta und Rückdrehung der Scheiben berechnen zu können
        
        % Berechnung des Grund-Punktes fuer das Segment der 3 Tendons
        for i=1:3
            w = (i*2*pi)/number_of_tendons; % Winkelschrittgroesse aus Tendonanzahl
            p_tendon(:,i) = [TACR.diskPitchRadius(j)*cos(w); TACR.diskPitchRadius(j)*sin(w); 0]; % Punkt der Tendons für z=0 aus Winkel(w) berechnen
        end
        delta = TACR.diskPitchRadius(j) * cos(phi_total); % delta in Abhaengigkeit zur gesamten Verdrehung um die z-Achse
        % Rotationsmatrix fuer Drehung um phi ist für Segment fest
        Rot_z = [cos(phi)    -sin(phi)    0;
                 sin(phi)     cos(phi)    0;
                    0           0         1];
        % Da beim ersten Segment der Verbindungsraum zum Vorgangssegment
        % wegfaellt wird die Schleife der Positionsberechnung einmal weniger
        % durchlaufen. 
        if j == 1
            s_num = TACR.ndisks(j);
        else
            s_num = TACR.ndisks(j)+1;
        end
        % Schleife ueber die jeweils auftretenden Abstaende(bei j>1 um eins großer)
        for s = 1:s_num
            % Drahtlaenge zur aktuellen Scheibe (s-1, da erste Scheibe auf dem boden liegen soll)
            q_s = (q(j,1)/10)*(s-1);
            % Winkel der Scheibe
            theta = theta_old + q_s(1)/delta;
            % Radian für Scheibe berechnen, wenn j=1 sind die Abstaende, welche
            % aus Segementlaenge resultieren nicht auch auf die Verbindung
            % zum Vorgangssegment verteilt. Beim Endstueck wird die dicke
            % des Endstückes abgezogen um die Segmentlaenge zu erzeugen wenn
            % die letzte Scheibe verschoben wird. 
            switch j
                case 1
                    radian = ((TACR.segmentLength(j)-3)/(TACR.ndisks(j)-1))*(s-1);
                case number_of_segments
                    radian = ((TACR.segmentLength(j)-4.3)/(TACR.ndisks(j)))*(s-1);
                otherwise
                    radian = (TACR.segmentLength(j)/(TACR.ndisks(j)))*(s-1);
            end
            alpha = theta_old - theta;
            % Wenn die Werte fuer q1/2/3=0 wird alpha 0 und somit r zu Inf
            % bei den Position treten folgend nur noch NaN-Werte auf
            if abs(alpha) < 0.0000001 % double-Wert, also wird über Minimum überprüft
                diskPointTemp = [0; 0; radian]; % Im spezialfall wird die Bogenlaenge als z-Wert gesetzt, da der Roboter senkrecht verlauft
            else
                r = radian/alpha;% Radius
                % Abstand des Punkts in x,y Ebene zum Ursprung (z(die Hoehe) wird nur mit alpha und dem Radius berechnet)
                a = r - cos(alpha)*r;
                diskPointTemp = [cos(phi)*a; sin(phi)*a; sin(alpha)*r];
            end
            % Rotationsmatrix fuer Neigung in y wird berechnen
            Rot_y = [ cos(alpha)     0     sin(alpha);
                         0           1         0     ;
                     -sin(alpha)     0     cos(alpha)];
            Rot = rot_old*Rot_z*Rot_y; % Rotationsmatrix wird aus Endrotation des alten Segments und der neuen Rotation berechnet
            % Transformation fuer Rueckrotation um z aufstellen. 
            % Hierbei wird ein z-Einheitsvektor entsprechend der
            % Scheibendrehung rotiert. Mit diesem wird folgend eine
            % Transformationsmatrix um den Vektor um den gesamten Winkel
            % phi berechnet.
            VDreh = Rot*[0; 0; 1]; % z-Einheitsvektor
            VDreh = VDreh/norm(VDreh); % in Scheibenausrichtung rotieren und normieren. 
            Trans = [VDreh(1)^2*(1-cos(phi_total))+cos(phi_total)                  VDreh(1)*VDreh(2)*(1-cos(phi_total))+VDreh(3)*sin(phi_total)    VDreh(1)*VDreh(3)*(1-cos(phi_total))-VDreh(2)*sin(phi_total); 
                     VDreh(1)*VDreh(2)*(1-cos(phi_total))-VDreh(3)*sin(phi_total)  VDreh(2)^2*(1-cos(phi_total))+cos(phi_total)                    VDreh(2)*VDreh(3)*(1-cos(phi_total))+VDreh(1)*sin(phi_total); 
                     VDreh(1)*VDreh(3)*(1-cos(phi_total))+VDreh(2)*sin(phi_total)  VDreh(2)*VDreh(3)*(1-cos(phi_total))-VDreh(1)*sin(phi_total)    VDreh(3)^2*(1-cos(phi_total))+cos(phi_total)];
            
            Rot_Scheibe = Trans*Rot;% Die Transformation der Scheibe resultiert aus deren Rotation und der berechneten Ruecktransformation
            diskPointTemp = rot_old*diskPointTemp+points_old; % Die Punkte resuktieren aus dem Endpunkt des letzten Segements und den neuen rotierten Punkten

            if s==1 && j>1  % Sofern Scheiben nummer =1 und Segment groesser als eins ist wird nur die Verbindung zwischen den Segmenten berechnet und wird folglich nicht als Scheibe gespeichert
            else
                robotShape.diskRotation(speicher,1:9) = [Rot_Scheibe(1,:) Rot_Scheibe(2,:) Rot_Scheibe(3,:)];% Rotation Speichern
                % Da Endscheibe dicker, wird dies überprueft und wenn die
                % Endscheibe berechnet wird wird der Dickeunterschied von
                % 7.3-3 = 4.3 aufaddiert.
                if (j == number_of_segments && s == TACR.ndisks(j)+1) ||(j == 1 && s == TACR.ndisks(j) && number_of_segments==1 )
                    diskHeightEnd = 4.3;
                end
                % Die Punkte werden mit der Segmenthoehe, welche in
                % Scheibenrichtung rotiert wurde beaufschlagt und
                % gespeicher
                robotShape.diskPoints(speicher,1:3) = rot90(diskPointTemp+VDreh*(TACR.diskHeight+diskHeightEnd));
                % Die Punkte der Tendons werden aus den Scheibenkoordinaten
                % und den anfangs berechneten rotierten Tendonkoodinaten wiederrum mit der Scheibenhoehe berechnet
                robotShape.diskPoints(speicher,4:6) = rot90(Rot_Scheibe*p_tendon(:,1)+diskPointTemp+VDreh*(TACR.diskHeight+diskHeightEnd));
                robotShape.diskPoints(speicher,7:9) = rot90(Rot_Scheibe*p_tendon(:,2)+diskPointTemp+VDreh*(TACR.diskHeight+diskHeightEnd));
                robotShape.diskPoints(speicher,10:12) = rot90(Rot_Scheibe*p_tendon(:,3)+diskPointTemp+VDreh*(TACR.diskHeight+diskHeightEnd));
                speicher = speicher + 1;% Speicherstelle erhoehen
            end
        end
        % Enddaten des Segments werden fuer das naechstes Segment gespeichert
        points_old = diskPointTemp;
        rot_old= Rot;
        theta_old = theta;
    end
end
