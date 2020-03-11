#make signature subset
#OSKM, Human, Mouse, All, vnature reprogrammed
#gse102348 - complete course, mouse, OSKM
#gse114581 - incomplete course, mouse, OSKM
#gse127927 - 7F course (complete for 7F and YF), mouse
#gse28688 - incomplete, OSKM, human
#gse46321 - complete, OSKM, mouse
#gse59527 - incomplete, OSKM, mouse
#gse71255 - complete (t-test), OSKM, OEE, KEE, mouse
#gse81891 - incomplete (t-test), OSKM, OS(v)KM, human
#gse89455 - incomplete (28 day?), OSKM and various combinations, human


#cormat for signatures
signatures_list_main_BH = list(OSKM_signature[OSKM_signature$adj_pval < 0.05,],
                               full_signature[full_signature$adj_pval < 0.05,],
                               mouse_signature[mouse_signature$adj_pval < 0.05,],
                               human_ort_signature[human_ort_signature$adj_pval < 0.05,],
                               all_signature[all_signature$adj_pval < 0.05,])
names(signatures_list_main_BH) = c("OSKM", "full", "mouse", "human_ort", "all")
cormat_maker(gse_list_main)



