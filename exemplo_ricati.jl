using SpecialFunctions
using QuadGK
using Plots

#
# Calcula a discretização do termo não linear
#
# φ(t) 
#
# Para a equação de Ricatti
#
# y²(t) = φ(t) = y'(t) + 2/t² 
#  y(1) = 1
#
# que tem solução analítica y(t)=1/t
#
#

#
# O vetor de valores de y(t) é Y e a função φ(y) 
# é o termo não linear que opera em y(t)
#
function Calcula_a(Y::Vector,φ::Function)

    # Número de informações discretas em y
    ny = length(Y)

    # Aloca a0 e a1
    a0 = zeros(ny)
    a1 = zeros(ny)

    # Calcula os valores de a1
    for i=1:ny-1

        # Intervalo
        ΔY = Y[i+1] - Y[i]

        # Coeficiente a1
        a1[i] = (φ(Y[i+1])-φ(Y[i]))/ΔY

        # Integral de φ(y) no intervalo
        integral = quadgk(φ,Y[i],Y[i+1])[1]

        # Coeficiente a0
        a0[i] = integral/(ΔY) - a1[i]*(Y[i+1]^2 - Y[i]^2)/(2*ΔY)

        # Particularizado para φ(y) = y^2
        # a0nal = (Y[i+1]^3 - Y[i]^3)/(3*ΔY) - a1[i]*(Y[i+1]^2 - Y[i]^2)/(2*ΔY)

    end

    # Retorna os vetores de coeficientes
    return a0,a1

end


#
# Calcula os coeficientes c 
#
function Calcula_c(a0::Vector, a1::Vector)

   # Inicializa os vetores com os coeficientes c
   c0 = zero(a0)
   c1 = zero(a1)

   for i in LinearIndices(c0)

       c0[i] = a0[i] - sum(c0[1:i-1])
       c1[i] = a1[i] - sum(c1[1:i-1])

   end

   # Retorna os coeficientes
   return c0,c1
end


# Calcula a constante de integração c para o exemplo do Ricatti
function Calcula_c_l(a0l::T,a1l::T,t0::T,y0::T) where T


    # Cte -a1l*t0
    cte = -a1l*t0

    return  exp(cte)*(y0 + a0l/a1l -2/t0 - 2*a1l*exp(-cte)*expint(cte))

end

# Calcula a resposta y_t(t) em um intervalo
function Calcula_y_l(t::T,a0l::T,a1l::T,c_l::T) where T

    # Cte -a1l*t
    cte = -a1l*t

    -(a0l/a1l) + (2/t) + (2*a1l*exp(-cte)*expint(cte)) + c_t*exp(-cte)

end

#
# Função para fazer a visualização da função aproximada φ(y)
#
function Aproximacao(c0::Vector,c1::Vector,Y::Vector,y::Float64)

    # Inicializa a saída da função
    saida = 0.0

    # Bem idiota, mas vamos lá
    contador = 1
    for vy in Y

        # Verifica se y está abaixo do valor de discretização atual,
        # para ligarmos os H's
        if vy<y
            saida += c0[contador] + c1[contador]*y
            contador += 1
        else
            break
        end

    end

   return saida


end


function main()

  # Cria a função φ(y)
  φ(y) = y^2

  # Vamos fazer uma discretização para o Y
  Y = collect(0.0:1.0:10)

  # Calcula os coeficientes a
  a0,a1 = Calcula_a(Y,φ)

  # Calcula os coeficientes c
  c0,c1 = Calcula_c(a0,a1)

  # Monta uma aproximação para a função
  ap(y) = Aproximacao(c0,c1,Y,y)

  # Monta os gráficos
  display(plot(Y,φ.(Y),label="Original"))
  display(plot!(Y,ap.(Y),label="Aproximada"))
  
  φ.(Y) .- ap.(Y)


  # Agora vamos calcular a resposta para o nosso problema,
  # primeiro intervalo
  t_l = 1.0
  y_l = 1.0

  # Coeficiente de integração
  c_l = Calcula_c_l(a0[1],a1[1],t_l,y_l) 

  # Resposta nesse intervalo
  y(t) = Calcula_y_l(t,a0[1],a1[1],c_l)

end