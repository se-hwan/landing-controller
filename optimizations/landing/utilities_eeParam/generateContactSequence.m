function [contactSequence] = generateContactSequence(contactPhases)
    
    % INPUT:    contactPhases (n x 3) - n contact points, and specified number of flight and contact phases + the beginning order
    % OUTPUT:   cs (n x ~) - full contact sequence for all n legs

    n = numel(contactPhases(:,1));
    cs = cell(n, 1);
    for i = 1:n
        cs{i} = legContactSequence(contactPhases(i,:));
    end
    
    contactSequence = cs;

end

function [leg_cs] =  legContactSequence(contactPhase)
    
    contact = contactPhase(1); flight = contactPhase(2); order = contactPhase(3);

    % INPUT:    contact (1) - number of contact modes
    %           flight (1) - number of flight modes
    %           order (1 or 0) - specifies whether model begins in contact or flight first
    % OUTPUT:   cs (~), full contact sequence of leg
    
    if (abs(flight - contact) > 1)
        error("Impossible contact sequence, specify number of modes appropriately");
    end
        
    switch order
        case 1 % leg begins in contact
            if (contact >= flight)
                leg_cs = ones(1, contact + flight);
                for i = 1:numel(leg_cs)
                    if (mod(i,2) == 0)
                        leg_cs(i) = 0;
                    end
                end
            else
                error("Cannot condition leg to begin in contact with fewer phases than flight");
            end
        case 0 % leg begins in flight
            if (flight >= contact)
                leg_cs = zeros(1, contact + flight);
                for i = 1:numel(leg_cs)
                    if (mod(i,2) == 0)
                        leg_cs(i) = 1;
                    end
                end
            else
                error("Cannot condition leg to begin in flight with fewer phases than contact");
            end
        otherwise
            error("Invalid order, specify 1 for leg to begin in contact and 0 for leg to begin in flight");
    end
end
