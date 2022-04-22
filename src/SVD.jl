using LinearAlgebra: size, norm, SVD, I, svd!
import LinearAlgebra.svd
using BlockArrays
using Match

@enum ReducedSVD begin
	Full = 0
	Thin = 1
	Compact = 2
	Truncated = 3
end

# TODO use specialization to exploit dynamic dispatching
# TODO return SVD factorization
export svd
function svd(A::AbstractBlockMatrix{T}; k=Nothing, ϵ::Real=1f-9, version::ReducedSVD=Thin) where {T}
	_A = if size(A, 1) < size(A, 2)
			A'
		else
			A
		end

	U,V = @match version begin
		Thin => OneSidedBlockJacobi(_A, ϵ)
		Full => error("Not implemented yet")
		Compact => error("Not implemented yet")
		Truncated => error("Not implemented yet")
	end

	if k != Nothing
		# U[:,:] = sortslices(U, dims=2, rev=true, by=sum)
		∑ = vec(sum(U, dims=1) / sum(U))
		p = sortperm(∑, rev=true)
		U[:,:] = U[:,p]
		V[:,:] = V[p,:]
	end

	return U,V
end

export OneSidedBlockJacobi
function OneSidedBlockJacobi(A::AbstractBlockMatrix{T}, ϵ) where {T}
	V = BlockArray{T}(I, blocksizes(A, 2), blocksizes(A, 2))
	U = copy(A)
	N =  blocksize(A, 2)

	checks = true
	while checks
		checks = false

		for i in 1:N
			for j in i + 1:N
				Uᵢ = blocks(U)[:,i]
				Uⱼ = blocks(U)[:,j]

				J = find_rotation!(Uᵢ, Uⱼ, ϵ)

				if J != Nothing
					# U⁽ᵏ⁾ := U⁽ᵏ⁻¹⁾ J
					pivot!(Uᵢ, Uⱼ, J)

					# V⁽ᵏ⁾ := V⁽ᵏ⁻¹⁾ J
					Vᵢ = blocks(V)[:,i]
					Vⱼ = blocks(V)[:,j]
					pivot!(Vᵢ, Vⱼ, J)

					checks = true
				end
			end
		end
	end

	return U, V
end

function find_rotation!(Aᵢ, Aⱼ, ϵ)
	bᵢᵢ = Aᵢ' * Aᵢ
	bⱼⱼ = Aⱼ' * Aⱼ
	bᵢⱼ	= Aᵢ' * Aⱼ

	if norm(bᵢⱼ) >= ϵ * sqrt(sum(bᵢᵢ .* bⱼⱼ))
		B = [bᵢᵢ bᵢⱼ; bᵢⱼ' bⱼⱼ]
		J, _, _ = svd!(B)
		return J
	else
		return Nothing
	end
end

# TODO change name for pivoting?
function pivot!(Aᵢ, Aⱼ, J)
	Aᵢⱼ = mortar([Aᵢ Aⱼ])
	J = PseudoBlockMatrix(J, [size(J, 1)], blocksizes(Aᵢⱼ, 2))

	Aᵢⱼ[:,:] = Aᵢⱼ * J
end