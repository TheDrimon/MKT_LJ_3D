import pygame as pg, numpy as np


def init(Type, n = 100, v = 1, di = 0):
    global Ox,Oy,Oz, Ovx,Ovy,Ovz, lax,lay,laz
    global Sx,Sy,Sz###
    
    if Type == 1:
        k = int(n**(1/3))
        
        a, b, c = np.meshgrid(np.linspace(xm+(xM-xm)/k/2, xM-(xM-xm)/(k)/2, num = k), np.linspace(ym+(yM-ym)/k/2, yM-(yM-ym)/k/2, num = k),np.linspace(zm+zM/k/2, zM-zM/k/2, num = k)) 
        Ox = a.flatten()
        Oy = b.flatten()
        Oz = c.flatten()
        
        #print("abc",a,b,c)
        
        alp = np.random.random(k*k*k) * 2 * np.pi
        bet = np.random.random(k*k*k) * 2 * np.pi
        Ovx = np.cos(bet) * np.sin(alp) * v
        Ovy = np.cos(bet) * np.cos(alp) * v
        Ovz = np.sin(bet) * v
        
        ax = np.zeros(len(Ox), dtype = np.float64)
        ay = np.zeros(len(Ox), dtype = np.float64)
        az = np.zeros(len(Ox), dtype = np.float64)
        Sx = np.zeros(len(Ox), dtype = np.float64)
        Sy = np.zeros(len(Ox), dtype = np.float64)
        Sz = np.zeros(len(Ox), dtype = np.float64)
        
        print("Sx",Sx)
        
        '''
        for i in range(len(Ox)):
            Dx = (Ox[i] - Ox)
            Dy = (Oy[i] - Oy)
            Dz = (Oz[i] - Oz)
            
            Dx = (Dx + (xM-xm)/2)%(xM-xm) - (xM -xm)/2
            Dy = (Dy + (yM-ym)/2)%(yM-ym) - (yM -ym)/2
            Dz = (Dz + (zM-zm)/2)%(zM-zm) - (zM -zm)/2
            
            Dx[i] = 1e99
            
            neRR = (r1*r1) / ( Dx*Dx + Dy*Dy + Dz*Dz ) # r1^2 / RR            
            neRR[i] = 0
            ##r1^2/6 = 4E       E = 215/24
            F = (2*np.power(neRR, 7) - np.power(neRR, 4)) * dt
                
            ax[i] = np.sum(F * Dx)
            ay[i] = np.sum(F * Dy)
            az[i] = np.sum(F * Dz)'''
        Ovx -= np.mean(Ovx)
        Ovy -= np.mean(Ovy)
        Ovz -= np.mean(Ovz)
        
        print((np.mean(Ovx)**2+np.mean(Ovy)**2+np.mean(Ovz)**2)**.5)
            	
        lax,lay,laz = ax,ay,az
    
    elif Type == 2:
        
        Ox=np.float64([100,100+di])
        Oy=np.float64([500,650])
        Oz=np.float64([200,200])
        Ovx=np.float64([0,0])
        Ovy=np.float64([0,-v]) * 2**.5
        Ovz=np.float64([0,0])
        
        lax = np.zeros(2,dtype = np.float64)
        lay = np.zeros(2,dtype = np.float64)
        laz = np.zeros(2,dtype = np.float64)
        Sx,Sy,Sz = lax,lay,laz###

#def detect(neRR):
#       return np.sum(neRR > 1.33333)
'''
def colided(neRR,i):
            global Hist
            if Hist[i][0]:
                if np.max(neRR) < 1:
                    Hist[i] = [0,Ox[i],Oy[i]]## not oz
            else:
                if np.max(neRR) > 1:
                    lll = Hist[i]
                    Hist[i] = [1,Ox[i],Oy[i]]
                    if lll[1] != -1:
                        print((lll[1]-Ox[i])**2 + (lll[2]-Oy[i])**2)
                        return (lll[1]-Ox[i])**2 + (lll[2]-Oy[i])**2

def colidedV(neRR,i):
            global Hist
            if Hist[i][0]:
                if np.max(neRR) < 1.333:
                    V = Ovx[i]*Ovx[i] + Ovy[i]*Ovy[i] + Ovz[i]*Ovz[i]
                    Hist[i] = [0,V,f*dt]
            else:
                if np.max(neRR) > 1.333:
                    lll = Hist[i]
                    V = Ovx[i]*Ovx[i] + Ovy[i]*Ovy[i] + Ovz[i]*Ovz[i]
                    Hist[i] = [1,V,f*dt]
                    if lll[1] != -1:
                        return ((V+lll[1])/2)**.5 * (Hist[i][2] - lll[2])
'''

nerrRconst = ( 0.96 ) ** (-2)# 0.88. 0.96

def colidedR(neRR,i):
            global Hist,Tcount
            if Hist[i][0]:
                if np.max(neRR) < nerrRconst:
                    ldj = Hist[i][-1]
                    #Hist[i] = [0,Sx[i],Sy[i],Sz[i],ldj]
                    Hist[i][0] =0
            else:
                if np.max(neRR) > nerrRconst:
                    lll = Hist[i]
                    #print(lll)
                    Hist[i] = [1,Sx[i],Sy[i],Sz[i],f]
                    if lll[1] != -1:
                        Tcount += [(Hist[i][-1]-lll[-1])*dt]
                        return ( (Hist[i][1]-lll[1])**2 +(Hist[i][2]-lll[2])**2 + (Hist[i][3]-lll[3])**2 )**.5 


#------------------------------------------------------------------------------------------------------------------------------------------------------------------#
def G_calc(dt): # единственный коректный
        global Ox,Oy,Oz, Ovx,Ovy,Ovz, lax,lay,laz,  Sx,Sy,Sz
        
        count = []
        
        global П
        П = 0
        
        ax = np.zeros(len(Ox),dtype = np.float64)
        ay = np.zeros(len(Ox),dtype = np.float64)
        az = np.zeros(len(Ox),dtype = np.float64)
        
        dx,dy,dz = ( Ovx + lax/2) * dt, ( Ovy + lay/2) * dt, ( Ovz + laz/2) * dt
        
        Sx += dx#
        Sy += dy#
        Sz += dz#
        
        
        
        Ox+=dx
        Oy+=dy
        Oz+=dz
        
        Ox = (Ox - xm)%(xM-xm) + xm
        Oy = (Oy - ym)%(yM-ym) + ym
        Oz = (Oz - zm)%(zM-zm) + zm
        
        for i in range(len(Ox)):
            Dx = (Ox[i] - Ox)
            Dy = (Oy[i] - Oy)
            Dz = (Oz[i] - Oz)
            
            Dx = (Dx + (xM-xm)/2)%(xM-xm) - (xM -xm)/2
            Dy = (Dy + (yM-ym)/2)%(yM-ym) - (yM -ym)/2
            Dz = (Dz + (zM-zm)/2)%(zM-zm) - (zM -zm)/2
            
            Dx[i] = 1e99
            
            neRR = (r1*r1) / ( Dx*Dx + Dy*Dy + Dz*Dz ) # r1^2 / RR            
            neRR[i] = 0
            
            П += np.sum(np.power(neRR, 6) - np.power(neRR, 3)) * (r1*r1 / 12)
            ##r1^2/6 = 4E       E = 215/24
            F = (2*np.power(neRR, 7) - np.power(neRR, 4)) * dt
                
            ax[i] = np.sum(F * Dx)
            ay[i] = np.sum(F * Dy)
            az[i] = np.sum(F * Dz)
            
            
            Bb = colidedR(neRR,i)
            if Bb:
                count += [Bb]
                   
        Ovx += (ax + lax) / 2
        Ovy += (ay + lay) / 2
        Ovz += (az + laz) / 2
        
        lax,lay,laz = ax,ay,az
        
        return count
#------------------------------------------------------------------------------------------------------------------------------------------------------------------#            

def loop(rep):
    global f
    f=0
    t = [pg.time.get_ticks(),0,0,0,0]
    k = 1
    
    global П
    П,K,E = 0,0,0
    E0 = 0
    
    п = np.zeros(700,dtype = np.float64)
    к = -np.ones(700,dtype = np.float64)*1000
    
    M = np.zeros(6, dtype = np.int64)
    global Sx,Sy,Sz
    #global lSx,lSy,lSz
    lSx,lSy,lSz = np.zeros((3,len(Ox)), dtype = np.float64)
    
    #SL = np.zeros(140, dtype = np.float64)
    SL=[0]
    Ra = 1
    Lpos = []    
    
    global Hist
    Hist = [[0,-1,-1,-1,-1]]*len(Ox)
    ans = []
    global Tcount
    Tcount = []
    
    L = np.zeros(125, dtype = np.float64)
    
    Sgm = []
    
    while True:
        t0 = pg.time.get_ticks()
        screen.fill((0,0,0))
        
        
        t1=pg.time.get_ticks()
        for i in range(rep):#######               <<========<<<======<<<<
            ans += G_calc(dt)
            f += 1
        t2=pg.time.get_ticks()
        
        if f == 1200:
            ans = []      
            #Sgm = []
            #SL = []
       
        for i in [[xm,yM+N//2,xm,ym],[xM,ym,xm-N//2,ym],[xM,yM+N//2,xM,ym-N//2],[xM,yM,xm,yM]]:
        	pg.draw.line(screen,([50,50,50]),i[0:2],i[2:4],N)
        
        for i in range(len(Ox)):
        	pg.draw.circle(screen,(0,0,0),(Ox[i],Oy[i]),1+r1/2.2 - Oz[i]//50,1)
        	pg.draw.circle(screen,(0,100,100*i/len(Ox)),(Ox[i],Oy[i]),r1/2.2 - Oz[i]//50)
        	
        pg.draw.circle(screen,(100,10,0),(Ox[i],Oy[i]),r1/2.2 - Oz[0]//50)
        	  
        t3=pg.time.get_ticks()
        
        if k < 5:
        	t[k] = pg.time.get_ticks()
        	k += 1
        else:
        	k = 0
        	
        ###############################################################################
        screen.blit(font.render("N = "+str(len(Ox)) + " fps = " + str(int(4000/(t[k-1]-t[k-5]))), 1, (255,0,0)),(20, 20))
        screen.blit(font.render("t = "+str(pg.time.get_ticks())+ "  " + str(100*(t1-t0)//(t[k-1]-t[k-2])) + "  " +str(100*(t2-t1)//(t[k-1]-t[k-2]))+ "  " + str(100*(t3-t2)//(t[k-1]-t[k-2])), 1, (255,0,0)),(20, 100))
        
        vv = Ovx*Ovx + Ovy*Ovy + Ovz*Ovz
        K = np.sum(vv)/2
        
        if  f < 2 * rep:
            E0 = K+П
            K0 = K
            П0 = П
 
        screen.blit(font.render("П = " + str(round(П)), 1, (255,0,0)),(20+900, 220+1500))
        screen.blit(font.render("К = " + str(round(K)), 1, (255,0,0)),(20+900, 280+1500))
        screen.blit(font.render("E = " + str(round(K+П,5)), 1, (255,0,0)),(20+900, 340+1500))
        screen.blit(font.render("dnE = " + str(round(100*(K+П)/E0-100,5))+"%", 1, (255,0,0)),(20+900, 400+1500))
        
######################RR
        hp = 100 ##### <==
        if f%(rep*10) == 0:
            #SL = np.roll(SL,1,0)
            if f%(rep*10*hp) == 0:
                lSx = Sx+0
                lSy = Sy+0
                lSz = Sz+0
            SL+=[np.mean((Sx-lSx)**2 + (Sy-lSy)**2 + (Sz-lSz)**2)]
            #SL+=[np.log(np.mean(Sx*Sx + Sy*Sy + Sz*Sz))+ 0/np.log(dt*f)]
            
            
            
        if len(SL):
            dbbs = -500 / (max(SL)+1e-15)
            pg.draw.line(screen,([190,50,50]),(800,1000),(1500,1000),2)
            #pg.draw.line(screen,([130,150,50]),(800,1000),(1500,1000-500),2)
            #pg.draw.line(screen,([130,150,50]),(800,1000),(1500,1000-500),2)
            ssl = 0
            for i in range(len(SL)-1):
                Ii = 8
                #pg.draw.line(screen,([200,30,120]),(800+Ii*np.log(i)*700/len(SL),1000+(SL[i]-ssl)*dbbs),(800+Ii*np.log(1+i)*700/len(SL),1000+(SL[i+1]-ssl)*dbbs),4)
                if (i+1)%hp:
                    pg.draw.line(screen,([200 * i/ len(SL),30,120* i/ len(SL)]),(800+i%hp*700/hp,1000+(SL[i]-ssl)*dbbs),(800+(1+i)%hp*700/hp,1000+(SL[i+1]-ssl)*dbbs),4)
            
            RRR = np.mean((Sx)*Sx + (Sy)*Sy + (Sz)*Sz)
            screen.blit(font.render("R^2 = " + str(round(RRR)), 1, (255,0,0)),(20+800, 460+560))
            screen.blit(font.render("(R^2)/t = " + str(round(RRR/f/dt)), 1, (255,0,0)),(20+800, 460+620))
            #screen.blit(font.render("lmd = " + str(round(RRR / (f * dt * 2 * (K / len(Ox))**0.5),6  )), 1, (255,0,0)),(20, 520))
        
######################E
        п[0] =-200 -500 * (П - min(K0, П0))/(abs(K0)+abs(П0))
        к[0] =-200 -500 * (K - min(K0, П0))/(abs(K0)+abs(П0))
        
        pg.draw.line(screen,([190,50,50]),(100,2400),(1500,2400),2)
        pg.draw.line(screen,([190,50,50]),(800,1300),(1500,1300),2)
        for i in range(len(п)-1):
            pg.draw.line(screen,([20,30,120]),(100+2*i,2400+п[i]),(100+2*i,2400+п[i+1]),4)
            pg.draw.line(screen,([120,30,20]),(100+2*i,2400+к[i]),(100+2*i,2400+к[i+1]),4)
            pg.draw.line(screen,([120,60,120]),(100+2*i,2400+к[i]+п[i]),(100+2*i,2400+к[i+1]+п[i+1]),4)
            pg.draw.line(screen,([20,150,20]),(800+i,1300+.01*(к[i]+п[i]+1000)/dt**1.5),(800+i,1300+.01*(к[i+1]+п[i+1]+1000)/dt**1.5),4)
            
        п = np.roll(п,1,0)
        к = np.roll(к,1,0)
        
######################e^-x2
        if f == rep * 300:
            M = np.zeros(400, dtype = np.int64)
            #Sx = 0
            #Sy = 0
            #Sz = 0
            #SL = []
            Tcount = []
        
        for v2 in np.int32(Ovx*10):
            if 0 < v2 < len(M):
                M[v2] += 1
         
        for i in range(len(M)):
            bl = i**2/300
            if abs(np.log(M[i]/(f-rep*100)/len(Ox))+10) < 6:
                pg.draw.line(screen,([120,80,120]),(100+bl, 1100-(np.log(M[i]/(f-rep*100)/len(Ox))+10)*150000/2000+0/len(Ox)),(100+bl,1300),2)
         
######################L<D
        if f%(rep*10) == 0 and len(ans):
            L = np.roll(L,1,0)
            L[0] = -np.mean(ans)
        
        #pg.draw.line(screen,([130,150,50]),(800,400+(L[0])*dbbs),(1500,400),2)
        
        if len(L):
            dbbs = 200
            pg.draw.line(screen,([190,50,50]),(700,400),(1200,400),2)
            mmm = max(L)-min(L) + 1e-15
            
            for i in range(len(L)-1):
                pg.draw.line(screen,([20,130,120]),(700+i*4,400+(L[i]-L[-1])*dbbs/mmm),(700+i*4+4,400+(L[i+1]-L[-1])*dbbs/mmm),4)
        screen.blit(font.render("<l> = " + str(round(np.mean(ans),2)), 1, (255,0,0)),(20+700, 160+260))
           
        screen.blit(font.render("<T> = " + str(round(np.mean(Tcount),2)), 1, (255,0,0)),(20+700, 580-100))
        sVV = np.mean(vv**0.5)# <v>
        screen.blit(font.render("<v>  = " + str(round(sVV,4)), 1, (255,0,0)),(20, 560))
        screen.blit(font.render("<v'> = " + str(round((2*K/len(Ox))**.5/(3*np.pi/8)**.5 ,4)), 1, (255,0,0)),(20, 620))
        
        screen.blit(font.render("<l> = <v> <T> = " + str(round(sVV*np.mean(Tcount),2)), 1, (255,0,0)),(20+700, 580-30))
        screen.blit(font.render("<l> = 3D / <v> = " + str(round(RRR/(f*dt*sVV)/2,2)), 1, (255,0,0)),(20+800, 460+680))
        if len(Tcount):
            sgm = 5*RRR/f #1/(sVV/r1 * 2**.5 * np.mean(Tcount) * len(Ox)/1000)
            Sgm += [-sgm]
            
            pg.draw.line(screen,([80,50,80]),(20,400),(20+404,400),3)
            pg.draw.line(screen,([80,50,80]),(20,400-sgm),(20+404,400-sgm),3)
            for i in range(len(Sgm)-1):
                pg.draw.line(screen,([20,130,120]),(20+i*400/len(Sgm),400+(Sgm[i])),(20+(1+i)*400/len(Sgm)+4,400+(Sgm[i+1])),4)
        #print(Sgm)
            
            #screen.blit(font.render("sgm = " + str(round(sgm,10)), 1, (255,0,0)),(20, 640))
        #screen.blit(font.render("rmn = " + str(round((np.mean(Sgm)*np.pi)*.5/2,10)), 1, (255,0,0)),(20, 700))
        
        screen.blit(font.render("D | D'  = " + str(round(2*(1/(len(Ox)/1000*RRR/(f*dt*sVV)/r1/2)/np.pi)**0.5,5)), 1, (255,0,0)),(200, 1670))
        screen.blit(font.render('''T | D" = ''' + str(round(2*(1/(len(Ox)/1000*sVV*np.mean(Tcount)/r1)/np.pi)**0.5,5)), 1, (255,0,0)),(208, 1780))
        screen.blit(font.render('''l | D" = ''' + str(round(2*(1/(len(Ox)/1000*np.mean(ans)/r1)/np.pi)**0.5,5)), 1, (255,0,0)),(230, 1725))
        #screen.blit(font.render("D* = " + str(round(2**(1/6)*(1+r1*r1/Vo/144*1.8)**0.5,10)), 1, (255,0,0)),(200, 1900))
        screen.blit(font.render("f | Do* = " + str(round((2 / (1 + (1 + (2*(K+П) / len(Ox)) /(r1*r1/48))**0.5))**(1/6),5)), 1, (255,0,0)),(200, 1850))
        screen.blit(font.render("f | D*   = " + str(round(2**(2/3) / (1 + (1 + (2*(K+П) / len(Ox)) /(r1*r1/48))**0.5)**(1/6),5)), 1, (255,0,0)),(200, 1900))
        
 ###################### (x,y,z,t)[0] 
 
        pg.draw.circle(screen,(30,40,40),(1000, 150),10)
        pg.draw.circle(screen,(30,40,40),(1000, 150),100,2)
        
        Ra = max( Ra, ( Sx[0]**2 + Sy[0]**2 ) ** .5 )
        if f%rep*30 == 0:
            Lpos += [np.array((Sx[0],Sy[0]))]
            if len(Lpos) >= 1000:
                NLpos = []
                for i in range(500):
                    NLpos += [Lpos[i*2]]
                Lpos = NLpos
        
        for i in range(len(Lpos)-1):
            pg.draw.line(screen,([70,170,120]),(100*Lpos[i]/Ra + (1000,150)),(100*Lpos[i+1]/Ra+ (1000,150)),2)
        pg.draw.circle(screen,(30,40,40),(1000, 150),100*np.exp(round(np.log(Ra)))/Ra,1)
        
        pg.draw.circle(screen,(140,10,0),(1000+Sx[0]/Ra*100, 150+Sy[0]/Ra*100),r1/3)     
        
        #print("Sx",Sx[:10])
        #print("Sy",Sy[:10])
        #print("Sz",Sz[:10])
        
        ###############################################################################
        pg.display.update()


if __name__ == "__main__":
    pg.init()

    xm, xM = 1080,1580
    ym, yM = 20,520# высота
    zm, zM = 0, 500
    
    
    mon = np.int16([xM,yM])
    screen = pg.display.set_mode(mon)
    pg.display.set_caption('simulatiOpn')
    font = pg.font.SysFont(None, 96)
    N=20
    
    
    # П - r1^2 / 24 =9,375
    # K - 1/2
    global dt, Vo
    Vo = 15   #<=
    dt = np.float64(.01) # <=
    r1=30
    init(1,215 * 5.0001, Vo,0* r1 * 1.24)# KT / Epsilon = v*48/r1^2.         2.66.     N=60 <=> nu=p = 0,1 - устарело
    print((Vo*Vo*2/3) /(r1*r1/48))# KT / Epsilon ≈ 2.15
    loop(1)#r1/Epsilon^.5 
    # mv^2 / epsilon = 12 * v^2 / r1^2
    # epsilon = r1^2 /12 