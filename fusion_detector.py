def df_handler(df):  
# here is the function that check if there are intersections between the domainlist we send in and the domain result of these seqs.
# it returns the longest result between the rpsblast_result's qstart and qend
    if not df.empty:
        length = 0
        for index, row in df.iterrows():
            if row.values[12] in domain_list:
                if row.values[7] - row.values[6] > length:
                    length = row.values[7] - row.values[6]
                    title = row.values[12]
                    start = row.values[6]
                    end = row.values[7]
        try:
            return title, start, end
        except:
            return None
        
def panduan(a, b, c, d):   # here is the function that check the intersection of each scale. and also a < b, c < d.
    def neibu(a, b, c, d):
        if a >= c:
            return "inclusion"
        elif c >= b:
            return "not intersected",c - b  # if two domains are not intersected, the possibility that they are fusion proteins is very high
        else:
            if (b-c)/min(b-a, d-c) <= 0.1:  # if the scale of two domains' intersection only occupies less than 0.1, the possibility that they are fusion proteins is very high
                return "intersected",b-c, max(d-b, c-a)
    if d >= b:
        return neibu(a, b, c, d)
    else:
        return neibu(c, d, a, b)
    
def result_check(df):
    df_handler_res = df_handler(df)
    if df_handler_res == None:
        pass
    else:
        first_start, first_end = df_handler_res[1], df_handler_res[2]  # the domain scale of our interested domain
        x = 0
        for index, row in df.iterrows():
            start, end = row.values[6], row.values[7]
            result = panduan(first_start, first_end, start, end)
            if result == 'inclusion' or not result:
                pass
            else:
                # to print the info of the two domains
                if x == 0:
                    print(item.lstrip('.'))
                print(result)
                print(first_start, first_end, start, end)
                print(row.values)
                print('\n')
                x += 1

def main(enzyme, domain_list, item):  # domain_list: list_that_you_send_in, item:the loc of the orf that have our interested domain
    for item in li:
        if os.path.exists(f'{item}/S_{enzyme}.fasta'):   # check that if the short orf seq exists.
            seq_ids = []   # the list is used to save the seq_ids that exists in the file contains the short orf seq 
                           # because the seq in long seq file may not exist in short seq file as well.
            with open(f'{item}/S_{enzyme}.fasta') as f:  # open the short orf seq file
                for line in f.readlines():
                    if line.startswith('>'):
                        seq_id = line.strip('>').split('\t')[0]
                        seq_ids.append(seq_id)

            rps_file = item.replace('domain_result', 'rps') + '/L_UroS' # the loc of rps_file, it's my loc (rps replace domain_result)
            if not os.path.getsize(rps_file) == 0:
                df = pd.read_csv(f'{rps_file}',sep='\t',header=None)
                df[12] = df[12].apply(lambda x: x.split(',')[0])
                whole_df = [item for i in df.groupby(0) for item in i if isinstance(item,pd.core.frame.DataFrame)]
                for eve_df in whole_df:       # as there may be several seqs in a rps_file, so they must be analized seperatedly
                    if eve_df.iloc[0, 0] not in seq_ids:
                        continue
                    eve_df = eve_df[eve_df[12].str.contains('PRK')|eve_df[12].str.contains('COG')|eve_df[12].str.contains('pfam')|eve_df[12].str.contains('TIGR')]  
                    # PRK, COG, pfam, TIGR are the trusted sources by me.
                    result_check(eve_df)

if __name__ == '__main__':
    
    # the job of this script is to detect the domain fusion of our intersted domains and some other domains
    main(enzyme, domain_list, item)