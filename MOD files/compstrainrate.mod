  n)  k   k820309    �          2021.10.0   ��g                                                                                                          
       strain.f90 COMPSTRAINRATE                                                             
                                                                   
                          @                                        '�                    #COORD    #VEL    #ACC    #PRESS 	   #MASS 
   #DENS    #MAT    #TXX    #TXY    #TYY    #HXX    #HXY    #HYY    #ID                 �                                                                   #VECTOR_3D                       @                                        '                    #X                 � $                                                                     
  p          p            p                                       �                                                                  #VECTOR_3D                 �                                                    0              #VECTOR_3D                 �                                        	     H          
                �                                        
     P          
                �                                             X          
                �                                             	       `                 
  p          p          p            p          p                                       �                                             �          
                �                                             �       	   
                �                                             �       
   
                �                                             �          
                �                                             �          
                �                                             �          
                �                                             �                                @                                        '                    #X                 � $                                                                     
  p          p            p                          #         @                                                                #KERNAL_TYPE    #RIJ    #RIJ_X    #RIJ_Y    #H    #WIJ    #DWDX    #DWDY              
                                                                
  @                                             
                
                                                
                
                                                
                
                                                
                                                               
                                                                
                                                                
       &         @   @                                                               #A     #B !   #VECTOR_3D              
                                                               #VECTOR_3D              
                                           !                   #VECTOR_3D    &         @   @                                    "                           #A #   #B $   #VECTOR_2D              
                                           #                   #VECTOR_2D              
                                           $                   #VECTOR_2D    &         @   @                                    %                           #A &   #B '   #VECTOR_3D              
                                           &                   #VECTOR_3D              
                                           '                   #VECTOR_3D    &         @   @                                    (                           #A )   #B *   #VECTOR_2D              
                                           )                   #VECTOR_2D              
                                           *                   #VECTOR_2D    &         @   @                                    +                           #A ,   #B -   #VECTOR_2D              
                                           ,                   #VECTOR_2D              
                                           -                   
    p          p            p                          &         @   @                                    .                           #A /   #B 0   #VECTOR_3D              
                                           /                   #VECTOR_3D              
                                           0     
      &         @   @                                    1                           #A 2   #B 3   #VECTOR_2D              
                                           2                   #VECTOR_2D              
                                           3     
      &         @   @                                    4                           #A 5   #B 6   #VECTOR_3D              
                                           5                   #VECTOR_3D              
                                           6                   #VECTOR_3D    %         @   @                                    7                    
       #A 8   #B 9             
                                           8                   #VECTOR_3D              
                                           9                   #VECTOR_3D    %         @   @                                    :                    
       #A ;   #B <             
                                           ;                   #VECTOR_2D              
                                           <                   #VECTOR_2D    #         @     @                                    =                    #A >   #B ?                                                       >                   
     p          p            p                                    
                                           ?                   #VECTOR_3D    #         @     @                                    @                    #A A   #B B                                                       A                    #VECTOR_3D              
                                           B     
      #         @     @                                    C                    #A D   #B E                                                       D                    #VECTOR_3D              
                                           E                   #VECTOR_3D    #         @     @                                    F                    #A G   #B H                                                       G                   
     p          p            p                                    
                                           H                   #VECTOR_2D    #         @     @                                    I                    #A J   #B K                                                       J                    #VECTOR_2D              
                                           K                   #VECTOR_2D    %         @   @                                    L                    
       #A M             
  @                                        M     
      &         @   @                                    N                           #A O   #B P   #VECTOR_3D              
                                           O                   #VECTOR_3D              
                                           P     
      &         @   @                                    Q                           #A R   #B S   #VECTOR_2D              
                                           R                   #VECTOR_2D              
                                           S     
      #         @                                            T                 
   #KERNAL_TYPE U   #VISCO_TYPE V   #P W   #NTOTAL X   #H Y   #DP Z   #XMIN [   #YMIN \   #LENGTH ]   #HEIGHT ^             
  @                                        U                     
                                           V                     
D                                          W            �                       &                                           #PARTICLE              
                                           X                     
  @                                        Y     
                
                                           Z     
                
  @                                        [     
                
  @                                        \     
                
                                           ]     
                
                                           ^     
      #         @                                            _                    #KERNAL_TYPE `   #VISCO_TYPE a   #P b   #NTOTAL c   #H d   #DP e   #EP_VEL f   #XMIN g   #YMIN h   #LENGTH i   #HEIGHT j             
  @                                        `                     
                                           a                     
D                                          b            �       
                &                                           #PARTICLE              
                                           c                     
  @                                        d     
                
                                           e     
                
                                           f                                  &                                           #VECTOR_2D              
  @                                        g     
                
  @                                        h     
                
                                           i     
                
                                           j     
         �   "      fn#fn    �   H   J   PARTICLES    
  H   J   KERNAL #   R  �       PARTICLE+PARTICLES )   -  g   a   PARTICLE%COORD+PARTICLES "   �  _       VECTOR_3D+VECTORS $   �  �   a   VECTOR_3D%X+VECTORS '   �  g   a   PARTICLE%VEL+PARTICLES '   �  g   a   PARTICLE%ACC+PARTICLES )   e  P   a   PARTICLE%PRESS+PARTICLES (   �  P   a   PARTICLE%MASS+PARTICLES (     P   a   PARTICLE%DENS+PARTICLES '   U  �   a   PARTICLE%MAT+PARTICLES '     P   a   PARTICLE%TXX+PARTICLES '   i  P   a   PARTICLE%TXY+PARTICLES '   �  P   a   PARTICLE%TYY+PARTICLES '   	  P   a   PARTICLE%HXX+PARTICLES '   Y  P   a   PARTICLE%HXY+PARTICLES '   �  P   a   PARTICLE%HYY+PARTICLES &   �  P   a   PARTICLE%ID+PARTICLES "   I  _       VECTOR_2D+VECTORS $   �  �   a   VECTOR_2D%X+VECTORS    L	  �       COMKER+KERNAL *   �	  H   a   COMKER%KERNAL_TYPE+KERNAL "   8
  H   a   COMKER%RIJ+KERNAL $   �
  H   a   COMKER%RIJ_X+KERNAL $   �
  H   a   COMKER%RIJ_Y+KERNAL       H   a   COMKER%H+KERNAL "   X  H   a   COMKER%WIJ+KERNAL #   �  H   a   COMKER%DWDX+KERNAL #   �  H   a   COMKER%DWDY+KERNAL (   0  u      VECTOR_3D_MINUS+VECTORS *   �  _   a   VECTOR_3D_MINUS%A+VECTORS *     _   a   VECTOR_3D_MINUS%B+VECTORS (   c  u      VECTOR_2D_MINUS+VECTORS *   �  _   a   VECTOR_2D_MINUS%A+VECTORS *   7  _   a   VECTOR_2D_MINUS%B+VECTORS '   �  u      VECTOR_3D_PLUS+VECTORS )     _   a   VECTOR_3D_PLUS%A+VECTORS )   j  _   a   VECTOR_3D_PLUS%B+VECTORS '   �  u      VECTOR_2D_PLUS+VECTORS )   >  _   a   VECTOR_2D_PLUS%A+VECTORS )   �  _   a   VECTOR_2D_PLUS%B+VECTORS ;   �  u      VECTOR_PLUS_VECTOR_MULTIPLY_SCALAR+VECTORS =   q  _   a   VECTOR_PLUS_VECTOR_MULTIPLY_SCALAR%A+VECTORS =   �  �   a   VECTOR_PLUS_VECTOR_MULTIPLY_SCALAR%B+VECTORS 5   l  u      VECTOR_3D_MULTIPLY_SCALAR_3D+VECTORS 7   �  _   a   VECTOR_3D_MULTIPLY_SCALAR_3D%A+VECTORS 7   @  H   a   VECTOR_3D_MULTIPLY_SCALAR_3D%B+VECTORS 5   �  u      VECTOR_2D_MULTIPLY_SCALAR_2D+VECTORS 7   �  _   a   VECTOR_2D_MULTIPLY_SCALAR_2D%A+VECTORS 7   \  H   a   VECTOR_2D_MULTIPLY_SCALAR_2D%B+VECTORS &   �  u      CROSS_PRODUCT+VECTORS (     _   a   CROSS_PRODUCT%A+VECTORS (   x  _   a   CROSS_PRODUCT%B+VECTORS '   �  f      DOT_PRODUCT_3D+VECTORS )   =  _   a   DOT_PRODUCT_3D%A+VECTORS )   �  _   a   DOT_PRODUCT_3D%B+VECTORS '   �  f      DOT_PRODUCT_2D+VECTORS )   a  _   a   DOT_PRODUCT_2D%A+VECTORS )   �  _   a   DOT_PRODUCT_2D%B+VECTORS /     ^      VECTOR_3D_TO_SCALAR_3D+VECTORS 1   }  �   a   VECTOR_3D_TO_SCALAR_3D%A+VECTORS 1     _   a   VECTOR_3D_TO_SCALAR_3D%B+VECTORS )   x  ^      SCALAR_TO_VECTOR+VECTORS +   �  _   a   SCALAR_TO_VECTOR%A+VECTORS +   5  H   a   SCALAR_TO_VECTOR%B+VECTORS /   }  ^      VECTOR_3D_TO_VECTOR_3D+VECTORS 1   �  _   a   VECTOR_3D_TO_VECTOR_3D%A+VECTORS 1   :  _   a   VECTOR_3D_TO_VECTOR_3D%B+VECTORS /   �  ^      VECTOR_2D_TO_SCALAR_2D+VECTORS 1   �  �   a   VECTOR_2D_TO_SCALAR_2D%A+VECTORS 1   �  _   a   VECTOR_2D_TO_SCALAR_2D%B+VECTORS /   �  ^      VECTOR_2D_TO_VECTOR_2D+VECTORS 1   P  _   a   VECTOR_2D_TO_VECTOR_2D%A+VECTORS 1   �  _   a   VECTOR_2D_TO_VECTOR_2D%B+VECTORS $     _      SCALAR_ROOT+VECTORS &   m  H   a   SCALAR_ROOT%A+VECTORS .   �  u      DIVISION_BY_SCALAR_3D+VECTORS 0   *  _   a   DIVISION_BY_SCALAR_3D%A+VECTORS 0   �  H   a   DIVISION_BY_SCALAR_3D%B+VECTORS .   �  u      DIVISION_BY_SCALAR_2D+VECTORS 0   F   _   a   DIVISION_BY_SCALAR_2D%A+VECTORS 0   �   H   a   DIVISION_BY_SCALAR_2D%B+VECTORS !   �   �       STRAIN_RATE_SLIP -   �!  H   a   STRAIN_RATE_SLIP%KERNAL_TYPE ,   �!  H   a   STRAIN_RATE_SLIP%VISCO_TYPE #   <"  �   a   STRAIN_RATE_SLIP%P (   �"  H   a   STRAIN_RATE_SLIP%NTOTAL #   &#  H   a   STRAIN_RATE_SLIP%H $   n#  H   a   STRAIN_RATE_SLIP%DP &   �#  H   a   STRAIN_RATE_SLIP%XMIN &   �#  H   a   STRAIN_RATE_SLIP%YMIN (   F$  H   a   STRAIN_RATE_SLIP%LENGTH (   �$  H   a   STRAIN_RATE_SLIP%HEIGHT #   �$  �       STRAIN_RATE_NOSLIP /   �%  H   a   STRAIN_RATE_NOSLIP%KERNAL_TYPE .   �%  H   a   STRAIN_RATE_NOSLIP%VISCO_TYPE %   1&  �   a   STRAIN_RATE_NOSLIP%P *   �&  H   a   STRAIN_RATE_NOSLIP%NTOTAL %   '  H   a   STRAIN_RATE_NOSLIP%H &   c'  H   a   STRAIN_RATE_NOSLIP%DP *   �'  �   a   STRAIN_RATE_NOSLIP%EP_VEL (   N(  H   a   STRAIN_RATE_NOSLIP%XMIN (   �(  H   a   STRAIN_RATE_NOSLIP%YMIN *   �(  H   a   STRAIN_RATE_NOSLIP%LENGTH *   &)  H   a   STRAIN_RATE_NOSLIP%HEIGHT 