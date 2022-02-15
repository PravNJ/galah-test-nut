import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

#Data Prep - same as dbscan.ipynb

df = pd.read_csv("tsne_run_all_20220207.csv", index_col=0)
df_sobject_ids = pd.read_csv("GDR3_sobject_ids_non_null_and_null.csv", index_col=0)
df['sobject_id'] = df_sobject_ids['sobject_id'].values
df_cotar_best_emission_candidates = pd.read_csv("/home/pravn/Dropbox/Masters/emission-mnras-cotar/best_emission_candidates.csv",header=None)
df_cotar_best_emission_candidates.columns = ["sobject_id"]
df_merged = pd.merge(df, df_cotar_best_emission_candidates, on = ["sobject_id"], how='outer', indicator = True)
df_cotar = pd.merge(df, df_cotar_best_emission_candidates, on = ["sobject_id"], how='inner')

fig,ax = plt.subplots(figsize=(50,50))
sc = plt.scatter(df_cotar["vis_x"],df_cotar["vis_y"], color = 'hotpink', alpha = 1, edgecolors = 'dimgray', s=20)

norm = plt.Normalize(1,4)
cmap = plt.cm.RdYlGn

plt.title('t-SNE Map - Cotar et al Candidates Only',fontsize=20)
plt.xlabel('vis_x',fontsize=14)
plt.ylabel('vis_y',fontsize=14)

annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

def update_annot(ind):

    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), 
                           " ".join([df_cotar["sobject_id"].to_numpy()[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)


plt.show()
