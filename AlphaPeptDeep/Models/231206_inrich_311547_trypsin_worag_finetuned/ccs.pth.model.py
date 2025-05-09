import torch
import peptdeep.model.building_block as building_block
from peptdeep.model.model_shop import *
class Model(torch.nn.Module):
    """LSTM model for CCS prediction"""
    def __init__(self,
        dropout=0.1
    ):
        super().__init__()
        
        self.dropout = torch.nn.Dropout(dropout)
        
        hidden = 256

        self.ccs_encoder = (
            model_base.Encoder_26AA_Mod_Charge_CNN_LSTM_AttnSum(
                hidden
            )
        )

        self.ccs_decoder = model_base.Decoder_Linear(
            hidden+1, 1
        )

    def forward(self, 
        aa_indices, 
        mod_x,
        charges,
    ):
        x = self.ccs_encoder(aa_indices, mod_x, charges)
        x = self.dropout(x)
        x = torch.cat((x, charges),1)
        return self.ccs_decoder(x).squeeze(1)
