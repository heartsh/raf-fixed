//////////////////////////////////////////////////////////////////////
// PCT.ipp
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// PCT::PCT()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
PCT<RealT>::PCT(std::vector<SparseMatrix<RealT> *> &base_pairing_posteriors,
                std::vector<SparseMatrix<RealT> *> &alignment_posteriors,
                const bool toggle_verbose,
                const int num_iterations) :
    base_pairing_posteriors(base_pairing_posteriors),
    alignment_posteriors(alignment_posteriors),
    m(int(base_pairing_posteriors.size())),
    toggle_verbose(toggle_verbose),
    num_iterations(num_iterations)
{
    Assert(m*m == int(alignment_posteriors.size()), "Dimension mismatch.");
}

//////////////////////////////////////////////////////////////////////
// PCT::Accumulate()
//
// Accumulate probabilistic consistency transformation through z.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void PCT<RealT>::Accumulate(std::vector<RealT> &res, int x, int y)
{
    const SparseMatrix<RealT> XY = (x < y ? *alignment_posteriors[x*m+y] : SparseMatrix<RealT>(*alignment_posteriors[y*m+x], SparseMatrix<RealT>::TRANSPOSE));
    const SparseMatrix<RealT> YX = (y < x ? *alignment_posteriors[y*m+x] : SparseMatrix<RealT>(*alignment_posteriors[x*m+y], SparseMatrix<RealT>::TRANSPOSE));
    const SparseMatrix<RealT> &YY = *base_pairing_posteriors[y];
    const SparseMatrix<RealT> YYT(YY, SparseMatrix<RealT>::TRANSPOSE);
    
    Assert(XY.GetNumCols() == YY.GetNumRows(), "Dimension mismatch.");
    Assert(YY.GetNumCols() == YX.GetNumRows(), "Dimension mismatch.");
    Assert(XY.GetNumRows() == YX.GetNumCols(), "Dimension mismatch.");
    Assert(XY.GetNumRows() * YX.GetNumCols() == int(res.size()), "Dimension mismatch.");

    const int LX = XY.GetNumRows()-1;
    const int LY = XY.GetNumCols()-1;

    /*
    for (int i = 1; i <= LX; i++)
    {
        for (const SparseMatrixEntry<RealT> *XYiter = XY.GetRowBegin(i); XYiter != XY.GetRowEnd(i); ++XYiter)
        {
            const int k = XYiter->column;
            for (const SparseMatrixEntry<RealT> *YYiter = YY.GetRowBegin(k); YYiter != YY.GetRowEnd(k); ++YYiter)
            {
                const RealT XYxYY = XYiter->value * YYiter->value;
                const int l = YYiter->column;
                for (const SparseMatrixEntry<RealT> *YXiter = YX.GetRowBegin(l); YXiter != YX.GetRowEnd(l); ++YXiter)
                {
                    const int j = YXiter->column;
                    res[i * (LX+1) + j] = std::max(res[i * (LX+1) + j], XYxYY * YXiter->value);
                }
            }
            for (const SparseMatrixEntry<RealT> *YYiter = YYT.GetRowBegin(k); YYiter != YYT.GetRowEnd(k); ++YYiter)
            {
                const RealT XYxYY = XYiter->value * YYiter->value;
                const int l = YYiter->column;
                for (const SparseMatrixEntry<RealT> *YXiter = YX.GetRowBegin(l); YXiter != YX.GetRowEnd(l); ++YXiter)
                {
                    const int j = YXiter->column;
                    res[i * (LX+1) + j] = std::max(res[i * (LX+1) + j], XYxYY * YXiter->value);
                }
            }
        }
    }
    */
    
    for (int i = 1; i <= LX; i++)
    {
        for (const SparseMatrixEntry<RealT> *XYiter = XY.GetRowBegin(i); XYiter != XY.GetRowEnd(i); ++XYiter)
        {
            const int k = XYiter->column;
            for (const SparseMatrixEntry<RealT> *YYiter = YY.GetRowBegin(k); YYiter != YY.GetRowEnd(k); ++YYiter)
            {
                const RealT XYxYY = XYiter->value * YYiter->value;
                const int l = YYiter->column;
                for (const SparseMatrixEntry<RealT> *YXiter = YX.GetRowBegin(l); YXiter != YX.GetRowEnd(l); ++YXiter)
                {
                    const int j = YXiter->column;
                    res[i * (LX+1) + j] += XYxYY * YXiter->value;
                }
            }
            for (const SparseMatrixEntry<RealT> *YYiter = YYT.GetRowBegin(k); YYiter != YYT.GetRowEnd(k); ++YYiter)
            {
                const RealT XYxYY = XYiter->value * YYiter->value;
                const int l = YYiter->column;
                for (const SparseMatrixEntry<RealT> *YXiter = YX.GetRowBegin(l); YXiter != YX.GetRowEnd(l); ++YXiter)
                {
                    const int j = YXiter->column;
                    res[i * (LX+1) + j] += XYxYY * YXiter->value;
                }
            }
        }
    }
    
}

//////////////////////////////////////////////////////////////////////
// PCT::Transform()
//
// Perform consistency transformation on posteriors.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void PCT<RealT>::Transform()
{
    for (int iter = 0; iter < num_iterations; iter++)
    {
        // compute new base pairing posteriors
        
        std::vector<SparseMatrix<RealT> *> new_base_pairing_posteriors(m, static_cast<SparseMatrix<RealT> *>(NULL));
        
        for (int x = 0; x < m; x++)
        {
            if (toggle_verbose)
            {
                WriteProgressMessage(SPrintF("Reestimating base-pairing posteriors: (%d)...", x+1));
            }

            /*
            std::vector<RealT> new_table = base_pairing_posteriors[x]->GetUnsparse();
            std::vector<RealT> old_mask = base_pairing_posteriors[x]->GetUnsparseMask();

            for (int y = 0; y < m; y++)
            {
                if (y == x) continue;
                Accumulate(new_table, x, y);
            }

            for (size_t i = 0; i < new_table.size(); i++)
            {
                if (!old_mask[i] && new_table[i] < RealT(0.1))
                    new_table[i] = RealT(0);
            }
            
            new_base_pairing_posteriors[x] = new SparseMatrix<RealT>(&new_table[0],
                                                                     base_pairing_posteriors[x]->GetNumRows(),
                                                                     base_pairing_posteriors[x]->GetNumCols(),
                                                                     RealT(0));
            */
            
            std::vector<RealT> new_table = base_pairing_posteriors[x]->GetUnsparse();
            
            for (int y = 0; y < m; y++)
            {
                if (y == x) continue;
                Accumulate(new_table, x, y);
            }

            new_table = new_table / RealT(m);

            /*
            for (size_t i = 0; i < new_table.size(); i++) if (new_table[i] < RealT(0.001)) new_table[i] = RealT(0);
            new_base_pairing_posteriors[x] = new SparseMatrix<RealT>(&new_table[0],
                                                                     base_pairing_posteriors[x]->GetNumRows(),
                                                                     base_pairing_posteriors[x]->GetNumCols(),
                                                                     RealT(0));

            */
            new_base_pairing_posteriors[x] = new SparseMatrix<RealT>(*base_pairing_posteriors[x], &new_table[0]);
        }

        // replace old posteriors

        for (size_t i = 0; i < base_pairing_posteriors.size(); i++)
        {
            delete base_pairing_posteriors[i];
            base_pairing_posteriors[i] = new_base_pairing_posteriors[i];
        }
    }

    if (toggle_verbose)
    {
        WriteProgressMessage("");
    }     
}
