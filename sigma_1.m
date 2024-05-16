function sig1 = sigma_1(z)
%   SIGMA1 Výpočet
    for i = 1:1:length(z)
        position = z(i,:);
        sig1(i,:) = position/sqrt(1 + norm(position));
    end
end

