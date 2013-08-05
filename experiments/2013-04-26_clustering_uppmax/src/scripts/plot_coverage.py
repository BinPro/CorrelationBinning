fig = plt.figure()
fig.suptitle("em and kmeans recall and precision for different number of samples")
ax1 = fig.add_subplot(121)
l1=ax1.plot([x[0] for x in recall_em], [x[1] for x in recall_em], "ro-", label="em_rec")
l2=ax1.plot([x[0] for x in recall_kmeans], [x[1] for x in recall_kmeans], "bx-",label="kmeans_rec")
ax1.set_ylabel("ratio")
ax1.set_xlabel("samples")
plt.ylim(0,1)
plt.xlim(0,50)
ax2 = fig.add_subplot(122,sharey=ax1,sharex=ax1)

l3=ax2.plot([x[0] for x in precision_em], [x[1] for x in precision_em], "ro--",label="em_prec")
l4=ax2.plot([x[0] for x in precision_kmeans], [x[1] for x in precision_kmeans], "bx--",label="kmeans_prec")
ax2.set_xlabel("samples")
ax1.legend(loc=4, ncol=1)
ax2.legend(loc=4, ncol=1)

fig.savefig("asdf.png")

