# function align_ylabels!(axs...)
#     max_extent = maximum(axs) do ax
#         tight_yticklabel_spacing!(ax)
#         ax.yticklabelspace[]
#     end
#     for ax in axs
#         ax.yticklabelspace = max_extent
#     end
#     return
# end

# f = Figure()
# ax1 = Axis(f[1, 1], ylabel = "Y Label")
# ax2 = Axis(f[2, 1], limits = (nothing, (1000, 9000)), ylabel = "Y Label")
# ax3 = Axis(f[3, 1], limits = (nothing, (0.0001, 0.0002)), ylabel = "Y Label")
# align_ylabels!(ax1, ax2, ax3)
# f