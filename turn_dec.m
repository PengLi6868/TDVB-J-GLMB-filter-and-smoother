function mk = turn_dec(m)

mk(1,:) = m(2,:).*sin(m(1,:));
mk(2,:) = m(2,:).*cos(m(1,:));