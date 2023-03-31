function u=fastSweep2DCase(u,f)

[m,n]=size(u);

for i =2:m-1
    for j =2:n-1
        
        A =min(u(i-1,j),u(i+1,j));
        B =min(u(i,j-1),u(i,j+1));
        all_ = sort([A,B]);%排序，第一个是小的
        a1 = all_(1);
        a2 = all_(2);
        
        xba = a1+f(i,j);
        if xba <= a2
            x = xba;
        else
            xba = 0.5*(a1+a2+sqrt(2*f(i,j)^2-(a1-a2)^2));
            if a2 <= xba
                x = xba;
            end
        end
        
        u(i,j) =min(x,u(i,j));
        
    end
end

for i =m-1:-1:2
    for j =n-1:-1:2
        
        A =min(u(i-1,j),u(i+1,j));
        B =min(u(i,j-1),u(i,j+1));
        all_ = sort([A,B]);
        a1 = all_(1);
        a2 = all_(2);
        
        xba = a1+f(i,j);
        if xba <= a2
            x = xba;
        else
            xba = 0.5*(a1+a2+sqrt(2*f(i,j)^2-(a1-a2)^2));
            if a2 <= xba
                x = xba;
            end
        end
        
        u(i,j) =min(x,u(i,j));
        
    end
end

for i =m-1:-1:2
    for j =2:n-1
        
        A =min(u(i-1,j),u(i+1,j));
        B =min(u(i,j-1),u(i,j+1));
        all_ = sort([A,B]);
        a1 = all_(1);
        a2 = all_(2);
        
        xba = a1+f(i,j);
        if xba <= a2
            x = xba;
        else
            xba = 0.5*(a1+a2+sqrt(2*f(i,j)^2-(a1-a2)^2));
            if a2 <= xba
                x = xba;
            end
        end
        
        u(i,j) =min(x,u(i,j));
        
    end
end

for i =2:m-1
    for j =n-1:-1:2
        
        A =min(u(i-1,j),u(i+1,j));
        B =min(u(i,j-1),u(i,j+1));
        all_ = sort([A,B]);
        a1 = all_(1);
        a2 = all_(2);
        
        xba = a1+f(i,j);
        if xba <= a2
            x = xba;
        else
            xba = 0.5*(a1+a2+sqrt(2*f(i,j)^2-(a1-a2)^2));
            if a2 <= xba
                x = xba;
            end
        end
        
        u(i,j) =min(x,u(i,j));
        
    end
end